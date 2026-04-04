#!/usr/bin/env julia

using JSON3
using StatsPlots
using Dates
using Plots

# --- Constantes do modelo (devem ser as mesmas usadas em run_sdp.jl) ---
const tmax = 30
const f_bins = 50

# --- Função para carregar e processar dados (agora em log10) ---
function load_and_process_data(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    script_dir = @__DIR__
    data_dir = joinpath(script_dir, "data")
    output_file = joinpath(data_dir, "SDP_$(fisio_type)_$(run_date).json")

    if !isfile(output_file)
        println("Erro: Arquivo de resultados '$output_file' para '$fisio_type' não encontrado.")
        println("Por favor, execute o script SDP.jl para '$fisio_type' e para a data '$run_date' primeiro.")
        return Dict{Float64, Float64}()
    end

    println("Carregando dados de: $output_file")
    data = JSON3.read(output_file, Dict)
    T_data = data["T"]

    proportions = Dict{Float64, Float64}()

    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents

    for mass in masses
        mass_str = string(round(mass; digits=2))

        if haskey(T_data, mass_str)
            hmax_actual = length(T_data[mass_str])
            T_matrix_for_mass = Array{Float64}(undef, hmax_actual, f_bins, tmax)

            for h_idx in 1:hmax_actual
                for f_idx in 1:f_bins
                    for t_idx in 1:tmax
                        T_matrix_for_mass[h_idx, f_idx, t_idx] = T_data[mass_str][h_idx][f_idx][t_idx]
                    end
                end
            end

            # --- Seleciona as últimas faixas de gordura (40–50) ---
            selected_rows = 40:50
            count_open_habitat = 0.0

            for i in 1:tmax
                for j in selected_rows
                    if T_matrix_for_mass[1, j, i] == 2.0
                        count_open_habitat += 1
                    end
                end
            end

            proportions[mass] = count_open_habitat / (tmax * length(selected_rows))
        else
            println("Aviso: Dados para a massa $mass_str não encontrados no arquivo JSON para $fisio_type.")
        end
    end

    return proportions
end

# --- Execução principal ---

if length(ARGS) < 7
    println("Uso: julia plot_sdp_fat.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed> <exp_start> <exp_end> <exp_step> <run_date>")
    println("Ex: julia plot_starvation.jl browsers grazers mixed 1 3 0.25 20251020")
    exit(1)
end

fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]
exp_start = parse(Float64, ARGS[4])
exp_end   = parse(Float64, ARGS[5])
exp_step  = parse(Float64, ARGS[6])
run_date  = ARGS[7]

# --- Carrega os dados ---
browser_open_data = load_and_process_data(fisio_type_browser, exp_start, exp_step, exp_end, run_date)
grazer_open_data  = load_and_process_data(fisio_type_grazer,  exp_start, exp_step, exp_end, run_date)
mixed_open_data   = load_and_process_data(fisio_type_mixed,   exp_start, exp_step, exp_end, run_date)

# --- Prepara os dados para o gráfico (eixo x em log10) ---
exponents = collect(exp_start:exp_step:exp_end)
masses_log10 = exponents

browser_proportions = [browser_open_data[10.0^x] for x in masses_log10]
grazer_proportions  = [grazer_open_data[10.0^x]  for x in masses_log10]
mixed_proportions   = [mixed_open_data[10.0^x]   for x in masses_log10]

# --- Plotagem ---
p = plot(masses_log10, browser_proportions,
    label = "Browser",
    xlabel = "Body Mass (log10) ",
    ylabel = "Proportion choosing open habitat",
    title = "Foraging Strategies per Body Mass (Fat Individuals)",
    lw = 2, marker = :circle, color = :blue,
    ylim = (0, 1),
    legend = :best,
    size = (800, 600))

plot!(p, masses_log10, grazer_proportions,
    label = "Grazer", lw = 2, marker = :diamond, color = :green)

plot!(p, masses_log10, mixed_proportions,
    label = "Mixed Feeder", lw = 2, marker = :star5, color = :orange)

# --- Salvar gráfico ---
script_dir = @__DIR__
results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

output_image_file = joinpath(results_dir, "foraging_by_mass_$(run_date)_fat.png")
savefig(p, output_image_file)
println("Gráfico salvo em $(output_image_file)")
