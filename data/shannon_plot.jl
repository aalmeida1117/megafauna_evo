#!/usr/bin/env julia

using JSON3
using Plots
using Statistics

# --- Constantes do modelo ---
const tmax = 30
const f_bins = 50

# --- Função para calcular o índice de Shannon ---
function calculate_shannon_index(proportions::Vector{Float64})
    epsilon = 1e-10
    props_adjusted = [max(p, epsilon) for p in proportions]
    total = sum(props_adjusted)
    props_normalized = props_adjusted ./ total
    shannon = -sum(p * log(p) for p in props_normalized)
    return shannon
end

# --- Função para carregar e processar dados (usando log10) ---
function load_and_process_data(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    script_dir = @__DIR__
    data_dir = joinpath(script_dir, "data")
    output_file = joinpath(data_dir, "SDP_$(fisio_type)_$(run_date).json")
    
    if !isfile(output_file)
        println("Erro: Arquivo '$output_file' não encontrado.")
        return Dict{Float64, Float64}()
    end
    
    println("Carregando: $output_file")
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
        end
    end
    
    return proportions
end

# --- Execução principal ---
if length(ARGS) < 7
    println("Uso: julia shannon_plot.jl <browser> <grazer> <mixed> <exp_start> <exp_end> <exp_step> <run_date>")
    println("Exemplo: julia shannon_plot.jl browsers grazers mixed 1 3 0.25 20250730")
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
browser_data = load_and_process_data(fisio_type_browser, exp_start, exp_step, exp_end, run_date)
grazer_data  = load_and_process_data(fisio_type_grazer,  exp_start, exp_step, exp_end, run_date)
mixed_data   = load_and_process_data(fisio_type_mixed,   exp_start, exp_step, exp_end, run_date)

# --- Prepara os dados ---
exponents = collect(exp_start:exp_step:exp_end)
masses_log10 = exponents  # eixo x será log10(massa)

browser_props = [browser_data[10.0^x] for x in masses_log10]
grazer_props  = [grazer_data[10.0^x]  for x in masses_log10]
mixed_props   = [mixed_data[10.0^x]   for x in masses_log10]

# --- Calcula o índice de Shannon ---
shannon_indices = Float64[]
for i in 1:length(masses_log10)
    props = [browser_props[i], grazer_props[i], mixed_props[i]]
    shannon = calculate_shannon_index(props)
    push!(shannon_indices, shannon)
end

# --- Cria o gráfico do índice de Shannon ---
max_shannon = log(3)

p = plot(masses_log10, shannon_indices, 
         label = "Shannon Index",
         xlabel = "log₁₀(Body mass)",
         ylabel = "Shannon Diversity Index",
         title = "Shannon Diversity Index per Body Mass",
         lw = 3, 
         marker = :circle, 
         markersize = 4,
         color = :purple,
         legend = :best,
         size = (900, 600))

# Linha do máximo teórico
hline!(p, [max_shannon], 
       label = "Maximum diversity (log(3))", 
       linestyle = :dash, 
       color = :red, 
       lw = 2)

# --- Salva o gráfico ---
script_dir = @__DIR__
results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

output_file = joinpath(results_dir, "shannon_index_$(run_date).png")
savefig(p, output_file)
println("\nGráfico salvo em: $(output_file)")

# --- Salva os dados em JSON ---
shannon_output = Dict(
    "log10_masses" => masses_log10,
    "shannon_indices" => shannon_indices,
    "run_date" => run_date
)

data_file = joinpath(results_dir, "shannon_data_$(run_date).json")
open(data_file, "w") do f
    JSON3.write(f, shannon_output)
end
println("Dados salvos em: $(data_file)")

# --- Estatísticas ---
println("\n--- Estatísticas ---")
println("Mínimo: $(round(minimum(shannon_indices), digits=4))")
println("Máximo: $(round(maximum(shannon_indices), digits=4))")
println("Média: $(round(mean(shannon_indices), digits=4))")
println("Máximo teórico: $(round(max_shannon, digits=4))")
