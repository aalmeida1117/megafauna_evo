#!/usr/bin/env julia
###############################################################
# analyze_SDP.jl
# Calcula e plota ΔZ/ΔM (comportamento) e ΔF/ΔM (fitness)
# Autor: Anna Almeida + ChatGPT + Claude
# Data: 2025-10-22
###############################################################

using JSON3
using Statistics
using Plots
using Dates

###############################################################
# ======== FUNÇÕES AUXILIARES ========
###############################################################

# --- Carrega arquivo JSON de saída do SDP ---
function load_SDP_results(fisio_type::String, run_date::String)
    script_dir = @__DIR__
    data_dir = joinpath(script_dir, "data")
    file = joinpath(data_dir, "SDP_$(fisio_type)_$(run_date).json")
    
    if !isfile(file)
        error("Arquivo não encontrado: $file")
    end
    
    println("→ Carregando resultados de: $file")
    return JSON3.read(file, Dict)
end

# --- Percentual de escolhas abertas (patch = 2) ---
function percent_open(Tmass)
    hmax = length(Tmass)
    f_bins = length(Tmass[1])
    tmax = length(Tmass[1][1])
    
    open_choices = sum(Tmass[2][y][t] == 2 for y in 1:f_bins, t in 1:tmax)
    total = f_bins * tmax
    
    return 100 * open_choices / total
end

# --- Fitness médio no último passo temporal ---
function mean_fitness(Smass)
    hmax = length(Smass)
    f_bins = length(Smass[1])
    tmax = length(Smass[1][1])
    
    return mean([Smass[h][y][tmax] for h in 1:hmax, y in 1:f_bins])
end

# --- Calcula deltas em log10 ---
function compute_deltas_log10(exponents, Z, F)
    masses = 10 .^ exponents
    ΔZΔM = diff(Z) ./ diff(masses)
    ΔFΔM = diff(F) ./ diff(masses)
    
    return (ΔZΔM, ΔFΔM)
end

# --- Processa um tipo fisiológico ---
function process_fisio_type(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    data = load_SDP_results(fisio_type, run_date)
    S = data["S"]
    T = data["T"]
    
    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents
    
    Z = Float64[]
    F = Float64[]
    
    for mass in masses
        mass_str = string(round(mass; digits=2))
        if haskey(T, mass_str) && haskey(S, mass_str)
            push!(Z, percent_open(T[mass_str]))
            push!(F, mean_fitness(S[mass_str]))
        else
            println("Aviso: massa $(mass_str) não encontrada no JSON de $(fisio_type).")
        end
    end
    
    ΔZΔM, ΔFΔM = compute_deltas_log10(exponents, Z, F)
    
    return (exponents, ΔZΔM, ΔFΔM)
end

###############################################################
# ======== EXECUÇÃO PRINCIPAL ========
###############################################################

if length(ARGS) < 7
    println("Uso: julia analyze_SDP.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed> <exp_start> <exp_end> <exp_step> <run_date>")
    println("Exemplo: julia analyze_SDP.jl browsers grazers mixed 1 3 0.25 20251022")
    exit(1)
end

fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]
exp_start = parse(Float64, ARGS[4])
exp_end   = parse(Float64, ARGS[5])
exp_step  = parse(Float64, ARGS[6])
run_date  = ARGS[7]

println("\n=== Processando dados ===")

# Processa os três tipos fisiológicos
exponents_b, ΔZΔM_browsers, ΔFΔM_browsers = process_fisio_type(fisio_type_browser, exp_start, exp_step, exp_end, run_date)
exponents_g, ΔZΔM_grazers, ΔFΔM_grazers   = process_fisio_type(fisio_type_grazer, exp_start, exp_step, exp_end, run_date)
exponents_m, ΔZΔM_mixed, ΔFΔM_mixed       = process_fisio_type(fisio_type_mixed, exp_start, exp_step, exp_end, run_date)

# Usa o eixo log10 (remove o primeiro ponto porque diff reduz o tamanho em 1)
masses_log10_b = exponents_b[2:end]
masses_log10_g = exponents_g[2:end]
masses_log10_m = exponents_m[2:end]

println("\n=== Gerando gráficos combinados ===")
println("Tamanhos: Browsers=$(length(ΔZΔM_browsers)), Grazers=$(length(ΔZΔM_grazers)), Mixed=$(length(ΔZΔM_mixed))")

# ============================================================
# GRÁFICO 1: ΔZ/ΔM (Comportamento) - TODAS AS CURVAS JUNTAS
# ============================================================
p1 = plot(masses_log10_b, ΔZΔM_browsers,
          label="Browsers",
          lw=2, 
          marker=:circle,
          color=:blue,
          xlabel="Mass (log10)",
          ylabel="ΔZ / ΔM (change in % open per kg)",
          title="Behavioral Sensitivity to Mass",
          legend=:best,
          size=(800, 600))

plot!(p1, masses_log10_g, ΔZΔM_grazers, 
      label="Grazers", 
      lw=2, 
      marker=:diamond,
      color=:green)

plot!(p1, masses_log10_m, ΔZΔM_mixed, 
      label="Mixed Feeders", 
      lw=2, 
      marker=:star5,
      color=:orange)

hline!(p1, [0], color=:black, linestyle=:dash, label="", alpha=0.5)

# ============================================================
# GRÁFICO 2: ΔF/ΔM (Fitness) - TODAS AS CURVAS JUNTAS
# ============================================================
p2 = plot(masses_log10_b, ΔFΔM_browsers,
          label="Browsers",
          lw=2,
          marker=:circle,
          color=:blue,
          xlabel="Mass (log10)",
          ylabel="ΔF / ΔM (change in fitness per kg)",
          title="Fitness Sensitivity to Mass",
          legend=:best,
          size=(800, 600))

plot!(p2, masses_log10_g, ΔFΔM_grazers, 
      label="Grazers", 
      lw=2,
      marker=:diamond,
      color=:green)

plot!(p2, masses_log10_m, ΔFΔM_mixed, 
      label="Mixed Feeders", 
      lw=2,
      marker=:star5,
      color=:orange)

hline!(p2, [0], color=:black, linestyle=:dash, label="", alpha=0.5)

# ============================================================
# Salvando os gráficos
# ============================================================
script_dir = @__DIR__
results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

output_behavior = joinpath(results_dir, "ΔZ_combined_$(run_date).png")
output_fitness = joinpath(results_dir, "ΔF_combined_$(run_date).png")

png(p1, output_behavior)
png(p2, output_fitness)

println("✓ Gráfico de comportamento (ΔZ) salvo: $(output_behavior)")
println("✓ Gráfico de fitness (ΔF) salvo: $(output_fitness)")
println("\n✅ Análises concluídas com sucesso para todos os fisiotipos.")