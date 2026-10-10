#!/usr/bin/env julia

using JSON3
using StatsPlots
using Dates
using Plots
#using LaTeXStrings   # <- adicionado para usar rótulos bonitos no eixo x

const tmax = 30
const f_bins = 50

# ============================================================
# Função para carregar e processar dados
# ============================================================
function load_and_process_data(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    script_dir = @__DIR__
    data_dir = joinpath(script_dir, "data")
    output_file = joinpath(data_dir, "SDP_$(fisio_type)_$(run_date).json")

    if !isfile(output_file)
        println("Erro: Arquivo '$output_file' não encontrado.")
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

            count_open_habitat = 0.0
            for i in 1:tmax
                for j in 1:f_bins
                    if T_matrix_for_mass[1, j, i] == 2.0
                        count_open_habitat += 1
                    end
                end
            end

            proportions[mass] = count_open_habitat / (tmax * f_bins)
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
    println("Uso: julia plot_sdp.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed> <exp_start> <exp_end> <exp_step> <run_date>")
    println("Exemplo: julia plot_sdp.jl browsers grazers mixed 0 4 0.1 20260310")
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
browser_open_data = load_and_process_data(fisio_type_browser, exp_start, exp_step, exp_end, run_date)
grazer_open_data  = load_and_process_data(fisio_type_grazer, exp_start, exp_step, exp_end, run_date)
mixed_open_data   = load_and_process_data(fisio_type_mixed, exp_start, exp_step, exp_end, run_date)

# ============================================================
# Imprime valores de massa e proporção para cada tipo
# ============================================================
println("\n========================================")
println("BROWSERS - Massa e Proporção de Open Habitat:")
println("========================================")
for (mass, prop) in sort(collect(browser_open_data))
    println("Massa: $(round(mass, digits=2)) kg | Proporção: $(round(prop, digits=4))")
end

println("\n========================================")
println("GRAZERS - Massa e Proporção de Open Habitat:")
println("========================================")
for (mass, prop) in sort(collect(grazer_open_data))
    println("Massa: $(round(mass, digits=2)) kg | Proporção: $(round(prop, digits=4))")
end

println("\n========================================")
println("MIXED - Massa e Proporção de Open Habitat:")
println("========================================")
for (mass, prop) in sort(collect(mixed_open_data))
    println("Massa: $(round(mass, digits=2)) kg | Proporção: $(round(prop, digits=4))")
end
println("========================================\n")

# ============================================================
# Prepara dados ordenados para plotagem (em log10)
# ============================================================
exponents = collect(exp_start:exp_step:exp_end)
masses_log10 = exponents  # eixo x será log10(massa)

browser_proportions = [browser_open_data[10.0^x] for x in masses_log10]
grazer_proportions  = [grazer_open_data[10.0^x]  for x in masses_log10]
mixed_proportions   = [mixed_open_data[10.0^x]   for x in masses_log10]

# ============================================================
# Configuração global de estilo (fonte, tamanhos, etc.)
# ============================================================
default(
    fontfamily = "Helvetica",
    guidefont = font(18),
    tickfont = font(15),
    legendfont = font(13),
    dpi = 600
)

# ============================================================
# Plotagem aprimorada
# ============================================================
browser_color = "#1B9E77"
grazer_color  = "#D95F02"
mixed_color   = "#7570B3"

# Define os ticks desejados (somente potências inteiras de 10)
xticks_values = [0, 1, 2, 3, 4]  # eixo está em log10(massa)
xticks_labels = ["10⁰", "10¹", "10²", "10³", "10⁴"]

p = plot(
    masses_log10, browser_proportions,
    xlabel = "Body mass (kg)",
    ylabel = "Proportion of choices for open habitat",
    color = browser_color,
    framestyle = :box,
    marker = :none,
    legend = :none,
    lw = 5,
    ylim = (0, 1),
    size = (800, 650),
    grid = false,
    xticks = (xticks_values, xticks_labels)   # ← só 10⁰–10⁴
)

plot!(p, masses_log10, grazer_proportions,
    marker = :none,
    lw = 5,
    color = grazer_color
)

plot!(p, masses_log10, mixed_proportions,
    marker = :none,
    lw = 5,
    color = mixed_color
)

#Linha horizontal em y=0.5
hline!(p, [0.5],
    lw =2,
    linestyle = :dash,
    color = :black
)
# ============================================================
# Salvando o gráfico
# ============================================================
script_dir = @__DIR__
results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

output_image_file = joinpath(results_dir, "foraging_by_mass_$(run_date).png")
savefig(p, output_image_file)
println("Gráfico salvo em $(output_image_file)")