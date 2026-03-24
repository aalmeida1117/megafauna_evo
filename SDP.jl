#!/usr/bin/env julia
using JSON3
using Dates
using LinearAlgebra
using ProgressMeter
using Base.Threads
# Carrega funções auxiliares usadas no cálculo de energia, volume gástrico etc.
include("allometric_functions.jl")
using .AllometricFunctions

# Função para limitar valor entre mínimo e máximo
function bc(val, min_val, max_val)
    return max(min(val, max_val), min_val)
end

# Função interpolate corrigida
function interpolate(yvalue, S_current_mass_array, patch, t)
    lowy = floor(Int, yvalue)
    highy = lowy + 1
    qy = yvalue - lowy

    lowy = max(1, min(lowy, size(S_current_mass_array, 2)))
    highy = max(1, min(highy, size(S_current_mass_array, 2)))

    S_low = S_current_mass_array[patch, lowy, t + 1]
    S_high = S_current_mass_array[patch, highy, t + 1]

    if lowy == highy
        return S_low
    else
        return qy * S_high + (1 - qy) * S_low
    end
end

# Função principal que roda a programação dinâmica estocástica para uma massa
function runSDP(mass, mass_dict)
    mass = round(mass; digits=2)
    max_mcal = 3
    tmax = 30
    #d = [0.05, 0.055]
    d = [prob_pred_closed(mass), prob_pred_open(mass)]
    gutrate = 0.6
    hmax = 2

    for h in 1:hmax
        prob_raw = mass_dict["$mass"]["$h"]["prob"]
        nbins = length(mass_dict["$mass"]["$h"]["gains"])
        mass_dict["$mass"]["$h"]["prob"] = reshape(prob_raw, nbins, nbins)
    end

    patchtree = Array{Int64}(undef, hmax, hmax)
    for i in 1:hmax
        patchtree[i, 1] = i
        patchtree[i, 2:hmax] = setdiff(collect(1:hmax), i)
    end

    numpathways = size(patchtree)[2]
    travelcosts = zeros(Float64, hmax, hmax) .+ AllometricFunctions.energy_cost(mass)
    travelcosts[diagind(travelcosts)] .= 0

    starving_bins = zeros(Int, hmax)
    f_bins = 50
    g_bins = 5
    ymax = 4870 * 0.02 * (mass^1.19)
    yc = 4870 * 0.1 * 0.02 * (mass^1.19)

    S = zeros(hmax, f_bins, tmax)
    T = zeros(hmax, f_bins, tmax)

    fat_bin_size = ymax / f_bins
    gut_bin_size = (AllometricFunctions.gut_volume_g(mass) * max_mcal) / g_bins
    a = 1 / (ymax - yc)
    b = -yc / (ymax - yc)

    for h in 1:hmax
        for y in 1:f_bins
            fat_energy = y * fat_bin_size
            S[h, y, tmax] = fat_energy > yc ? a * fat_energy + b : 0
        end
        starving_bins = sum(S[h, :, tmax] .== 0.0)
    end

    nbins = length(mass_dict["$mass"]["1"]["costs"])

    for t in (tmax-1):-1:1
        for h in 1:hmax
            for y in (starving_bins + 1):f_bins
                surv = zeros(numpathways)
                for k in 1:numpathways
                    patch = patchtree[h, k]
                    travel = travelcosts[h, patch]
                    st = 0.0
                    for g in 1:g_bins
                        for i in 1:nbins
                            for j in 1:nbins
                                fat_state = y * fat_bin_size
                                gut_state = g * gut_bin_size
                                gain = mass_dict["$mass"]["$patch"]["gains"][i]
                                loss = mass_dict["$mass"]["$patch"]["costs"][j]
                                prob = mass_dict["$mass"]["$patch"]["prob"][i, j]
                                gutpass = (gain + gut_state) * gutrate
                                yp = fat_state + gutpass - loss - travel
                                yp_bin = bc(yp / fat_bin_size, 1, f_bins)
                                Svalue = interpolate(yp_bin, S, patch, t)
                                st += prob * Svalue * (1 / g_bins)
                            end
                        end
                    end
                    surv[k] = max((1 - d[h]) * st, 0)
                end
                maxvalue = maximum(surv)
                max_indices = findall(x -> x == maxvalue, surv)
                strat = rand(max_indices)
                S[h, y, t] = maxvalue
                T[h, y, t] = strat
            end
        end
    end

    return S, T
end


# ============================================================
# Salvando resultados
# ============================================================
function convert_array_3d_to_nested_vector(A::Array{Float64, 3})
    return [[[A[i, j, k] for k in 1:size(A, 3)] for j in 1:size(A, 2)] for i in 1:size(A, 1)]
end

function prepare_output_dict(dict::Dict)
    out = Dict()
    for (k, arr) in dict
        if isa(arr, Array{Float64, 3})
            # MUDANÇA: Garante que a chave seja string com 2 casas decimais
            mass_str = string(round(k; digits=2))  # <- MODIFICADO
            out[mass_str] = convert_array_3d_to_nested_vector(arr)
        else
            @warn "O valor para a chave '$k' não é um Array{Float64, 3}. Tipo: $(typeof(arr)). Pulando."
        end
    end
    return out
end

# ============================================================
# Parâmetros via terminal (escala log10)
# ============================================================
if length(ARGS) < 4
    println("Uso: julia SDP.jl <fisio_type> <exp_start> <exp_end> <exp_step>")
    println("Exemplo: julia SDP.jl grazers 1 3 0.25")
    exit(1)
end

fisio_type = ARGS[1]
exp_start = parse(Float64, ARGS[2])
exp_end = parse(Float64, ARGS[3])
exp_step = parse(Float64, ARGS[4])

S = Dict()
T = Dict()
###
# Carrega dados de entrada
script_dir = @__DIR__
data_dir = joinpath(script_dir, "data")
#input_file = joinpath(data_dir, "$(fisio_type)_dict.json")
mu_files = filter(f -> startswith(f, "$(fisio_type)_mu") && endswith(f, "_dict.json"), readdir(data_dir))
# if !isfile(input_file)
#     println("Erro: Arquivo de entrada '$(input_file)' não encontrado.")
#     println("Por favor, certifique-se de que o arquivo JSON está em $(data_dir).")
#     exit(1)
# end

if isempty(mu_files)
    println("Erro: Nenhum arquivo mu encontrado em $(data_dir).")
    exit(1)
end
#mass_dict = JSON3.read(input_file, Dict)

# Cria lista de massas como potências de 10
exponents = collect(exp_start:exp_step:exp_end)
masses = 10 .^ exponents

println("\nMassas (10^x): ", exponents)
println("Valores numéricos de massa: ", round.(masses; digits=2))

# Cria barra de progresso
println("\nProcessando $(length(masses)) massas...")
#progress = Progress(length(masses), desc="Calculando SDP: ", barlen=50)


# MUDANÇA: Arredonda a massa antes de usar como chave
@threads for mu_file in mu_files
    # extrai o valor de mu do nome do arquivo
    
    mu_str = replace(mu_file, "$(fisio_type)_mu" => "", "_dict.json" => "")
    println("mu=$mu_str rodando na thread $(Threads.threadid())")

    # cada thread carrega e processa seu próprio arquivo
    mass_dict_local = JSON3.read(joinpath(data_dir, mu_file), Dict)

    S_local = Dict()
    T_local = Dict()
    
    # for mass in masses
    #     mass_rounded = round(mass; digits=2)  # <- NOVA LINHA
    #     S[mass_rounded], T[mass_rounded] = runSDP(mass, mass_dict)
    #     next!(progress, showvalues = [(:massa_atual, mass_rounded)])
    # end
    for mass in masses
        mass_rounded = round(mass; digits=2)
        S_local[mass_rounded], T_local[mass_rounded] = runSDP(mass, mass_dict_local)
        println("  [mu=$mu_str] massa $mass_rounded concluída (thread $(Threads.threadid()))")
    end

    #finish!(progress)

    # salva um arquivo de saída por mu
    output_file = joinpath(data_dir, "SDP_$(fisio_type)_mu$(mu_str).json")
    S_out = prepare_output_dict(S_local)
    T_out = prepare_output_dict(T_local)
    JSON3.write(output_file, Dict("S" => S_out, "T"=> T_out))
    println("✓ mu=$mu_str salvo em $(output_file)")

end


# println("\nPreparando dados para salvar...")
# S_out = prepare_output_dict(S)
# T_out = prepare_output_dict(T)

# # Salva saída
# output_file = joinpath(data_dir, "SDP_$(fisio_type)_$(Dates.format(now(), "yyyymmdd")).json")
# println("Salvando resultados...")
# JSON3.write(output_file, Dict("S" => S_out, "T" => T_out))
# println("\n✓ Programa terminado. Resultados salvos em $(output_file)")