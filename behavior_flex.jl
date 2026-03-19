#!/usr/bin/env julia

using JSON3
using StatsPlots
using Plots
using Dates

const tmax = 30
const f_bins = 50

# ============================================================
# Função para carregar e processar dados comparando reservas de gordura
# ============================================================
function load_and_process_data(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    script_dir = @__DIR__
    data_dir = joinpath(script_dir, "data")
    output_file = joinpath(data_dir, "SDP_$(fisio_type)_$(run_date).json")
    
    if !isfile(output_file)
        println("Erro: Arquivo '$output_file' não encontrado.")
        return Dict{Float64, Dict{String, Float64}}()
    end
    
    println("Carregando dados de: $output_file")
    data = JSON3.read(output_file, Dict)
    T_data = data["T"]
    
    proportions = Dict{Float64, Dict{String, Float64}}()
    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents
    
    for mass in masses
        mass_str = string(round(mass; digits=2))
        
        if haskey(T_data, mass_str)
            hmax_actual = length(T_data[mass_str])
            f_bins_actual = length(T_data[mass_str][1])
            tmax_actual = length(T_data[mass_str][1][1])
            
            T_matrix_for_mass = Array{Float64}(undef, hmax_actual, f_bins_actual, tmax_actual)
            
            for h_idx in 1:hmax_actual
                for f_idx in 1:f_bins_actual
                    for t_idx in 1:tmax_actual
                        T_matrix_for_mass[h_idx, f_idx, t_idx] = T_data[mass_str][h_idx][f_idx][t_idx]
                    end
                end
            end
            
            # Define faixas de reservas (fat bins)
            # Baixas reservas: bins 6 a 16
            # Altas reservas: bins 40 a 50
            low_reserve_bins = 6:16
            high_reserve_bins = 40:50
            
            
            # Verifica se os ranges são válidos
            if isempty(low_reserve_bins) || isempty(high_reserve_bins)
                println("Aviso: Massa $(mass_str) - ranges vazios. low=$(low_reserve_bins), high=$(high_reserve_bins)")
                continue
            end
            
            # Calcula proporção para BAIXAS reservas (fat bins 6-16)
            count_open_low = 0.0
            for h in 1:hmax_actual
                for f in low_reserve_bins
                    for t in 1:tmax_actual
                        if T_matrix_for_mass[h, f, t] == 2.0
                            count_open_low += 1
                        end
                    end
                end
            end
            prop_low = count_open_low / (length(low_reserve_bins) * tmax_actual)
            
            # Calcula proporção para ALTAS reservas (fat bins 40-50)
            count_open_high = 0.0
            for h in 1:hmax_actual
                for f in high_reserve_bins
                    for t in 1:tmax_actual
                        if T_matrix_for_mass[h, f, t] == 2.0
                            count_open_high += 1
                        end
                    end
                end
            end
            prop_high = count_open_high / (length(high_reserve_bins) * tmax_actual)
            
            # Diferença: altas reservas - baixas reservas
            diff = prop_high - prop_low
            
            proportions[mass] = Dict(
                "low" => prop_low,
                "high" => prop_high,
                "diff" => diff
            )
            
            println("Massa $(mass_str): baixa=$(round(prop_low, digits=3)), alta=$(round(prop_high, digits=3)), diff=$(round(diff, digits=3))")
        else
            println("Aviso: massa $(mass_str) não encontrada no JSON de $(fisio_type).")
        end
    end
    
    return proportions
end

# ============================================================
# Execução principal
# ============================================================
if length(ARGS) < 7
    println("Uso: julia plot_habitat_by_reserves.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed> <exp_start> <exp_end> <exp_step> <run_date>")
    println("Exemplo: julia plot_habitat_by_reserves.jl browsers grazers mixed 1 3 0.25 20251022")
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
browser_data = load_and_process_data(fisio_type_browser, exp_start, exp_step, exp_end, run_date)
grazer_data  = load_and_process_data(fisio_type_grazer, exp_start, exp_step, exp_end, run_date)
mixed_data   = load_and_process_data(fisio_type_mixed, exp_start, exp_step, exp_end, run_date)

# ============================================================
# Prepara dados ordenados para plotagem (em log10)
# ============================================================
exponents = collect(exp_start:exp_step:exp_end)
masses_log10 = exponents  # eixo x será log10(massa)

# Extrai as diferenças
browser_diff = [browser_data[10.0^x]["diff"] for x in masses_log10]
grazer_diff  = [grazer_data[10.0^x]["diff"]  for x in masses_log10]
mixed_diff   = [mixed_data[10.0^x]["diff"]   for x in masses_log10]

# ============================================================
# Plotagem
# ============================================================
println("\n=== Gerando gráfico ===")
p = plot(masses_log10, browser_diff, 
         label="Browser", 
         xlabel="Mass (log10)",
         ylabel="Difference in proportion (high - low reserves)",
         title="Habitat Selection: High vs Low Body Reserves",
         lw=2, marker=:circle, color=:blue,
         legend=:best, size=(800, 600))

plot!(p, masses_log10, grazer_diff, label="Grazer", lw=2, marker=:diamond, color=:green)
plot!(p, masses_log10, mixed_diff, label="Mixed Feeder", lw=2, marker=:star5, color=:orange)
hline!([0], color=:black, linestyle=:dash, label="", alpha=0.5)

# ============================================================
# Salvando o gráfico
# ============================================================
script_dir = @__DIR__
results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

output_image_file = joinpath(results_dir, "habitat_by_reserves_$(run_date).png")
png(p, output_image_file)

println("✓ Gráfico salvo em: $(output_image_file)")