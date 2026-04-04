#!/usr/bin/env julia

using JSON3
using Statistics
using Plots
using Dates

const tmax = 30
const f_bins = 50

# ============================================================
# Função para carregar dados e calcular médias de fitness
# ============================================================
function load_and_calculate_means(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    script_dir = @__DIR__
    data_dir = joinpath(script_dir, "data")
    output_file = joinpath(data_dir, "SDP_$(fisio_type)_$(run_date).json")
    
    if !isfile(output_file)
        println("Erro: Arquivo '$output_file' não encontrado.")
        return Dict{Float64, Dict{Int, Float64}}()
    end
    
    println("Carregando dados de: $output_file")
    data = JSON3.read(output_file, Dict)
    S_data = data["S"]
    
    means_dict = Dict{Float64, Dict{Int, Float64}}()
    
    # Cria lista de massas
    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents
    
    for mass in masses
        mass_str = string(round(mass; digits=2))
        
        if haskey(S_data, mass_str)
            means_dict[mass] = Dict{Int, Float64}()
            
            hmax_actual = length(S_data[mass_str])
            
            for patch in 1:hmax_actual
                # Reconstrói a matriz S para este patch
                f_bins_actual = length(S_data[mass_str][patch])
                tmax_actual = length(S_data[mass_str][patch][1])
                
                mat = zeros(Float64, f_bins_actual, tmax_actual)
                for f_idx in 1:f_bins_actual
                    for t_idx in 1:tmax_actual
                        mat[f_idx, t_idx] = S_data[mass_str][patch][f_idx][t_idx]
                    end
                end
                
                # Calcula a média de todos os valores
                means_dict[mass][patch] = mean(mat)
            end
        else
            println("Aviso: massa $(mass_str) não encontrada no JSON de $(fisio_type).")
        end
    end
    
    return means_dict
end

# ============================================================
# Execução principal
# ============================================================
if length(ARGS) < 7
    println("Uso: julia plot_mean_fitness.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed> <exp_start> <exp_end> <exp_step> <run_date>")
    println("Exemplo: julia plot_mean_fitness.jl browsers grazers mixed 1 3 0.25 20251022")
    exit(1)
end

fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]
exp_start = parse(Float64, ARGS[4])
exp_end   = parse(Float64, ARGS[5])
exp_step  = parse(Float64, ARGS[6])
run_date  = ARGS[7]

# ============================================================
# Carrega dados para os três tipos fisiológicos
# ============================================================
println("\n=== Carregando dados ===")
means_browsers = load_and_calculate_means(fisio_type_browser, exp_start, exp_step, exp_end, run_date)
means_grazers  = load_and_calculate_means(fisio_type_grazer, exp_start, exp_step, exp_end, run_date)
means_mixed    = load_and_calculate_means(fisio_type_mixed, exp_start, exp_step, exp_end, run_date)

# ============================================================
# Prepara dados ordenados para plotagem (em log10)
# ============================================================
exponents = collect(exp_start:exp_step:exp_end)
masses_log10 = exponents  # eixo x será log10(massa)

# Separa por patch (habitat)
grazers_patch1  = [means_grazers[10.0^x][1]  for x in masses_log10]
grazers_patch2  = [means_grazers[10.0^x][2]  for x in masses_log10]
browsers_patch1 = [means_browsers[10.0^x][1] for x in masses_log10]
browsers_patch2 = [means_browsers[10.0^x][2] for x in masses_log10]
mixed_patch1    = [means_mixed[10.0^x][1]    for x in masses_log10]
mixed_patch2    = [means_mixed[10.0^x][2]    for x in masses_log10]

# ============================================================
# Plotagem
# ============================================================
println("\n=== Gerando gráfico ===")
p = plot(masses_log10, grazers_patch1, 
         label="Grazers - Closed Habitat", 
         lw=2, 
         legend=:outerright, 
         size=(800, 600),
         xlabel="Mass (log10)",
         ylabel="Mean fitness",
         title="Mean Fitness by Body Mass and Habitat")

plot!(p, masses_log10, grazers_patch2, label="Grazers - Open Habitat", lw=2)
plot!(p, masses_log10, browsers_patch1, label="Browsers - Closed Habitat", lw=2, ls=:dash)
plot!(p, masses_log10, browsers_patch2, label="Browsers - Open Habitat", lw=2, ls=:dash)
plot!(p, masses_log10, mixed_patch1, label="Mixed - Closed Habitat", lw=2, ls=:dot)
plot!(p, masses_log10, mixed_patch2, label="Mixed - Open Habitat", lw=2, ls=:dot)

# ============================================================
# Salvando o gráfico
# ============================================================
script_dir = @__DIR__
results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

output_image_file = joinpath(results_dir, "mean_fitness_by_mass_$(run_date).png")
savefig(p, output_image_file)

println("✓ Gráfico salvo em: $(output_image_file)")