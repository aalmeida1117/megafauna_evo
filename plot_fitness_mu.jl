#!/usr/bin/env julia
"""
plot_fitness_mu.jl
------------------
Lê os arquivos SDP_<fisio_type>_mu<valor>.json gerados pelo SDP.jl e plota:
  1. Fitness máximo vs mu  (uma curva por tipo: grazer, browser, mixed feeder)
  2. Massa correspondente ao fitness máximo vs mu

Uso:
  julia plot_fitness_mu.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed>
  Exemplo: julia plot_fitness_mu.jl browsers grazers mixed
"""

using JSON3
using Statistics
using StatsPlots
using Plots

# ── Argumentos ────────────────────────────────────────────────────────────────

if length(ARGS) < 3
    println("Uso: julia plot_fitness_mu.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed>")
    println("Exemplo: julia plot_fitness_mu.jl browsers grazers mixed")
    exit(1)
end

fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]

script_dir = @__DIR__
data_dir   = joinpath(script_dir, "data")

# ── Funções auxiliares ────────────────────────────────────────────────────────

"""Extrai o valor numérico de mu do nome do arquivo."""
function parse_mu(filename, fisio_type)
    pattern = Regex("SDP_$(fisio_type)_mu(.+)\\.json")
    m = match(pattern, basename(filename))
    m === nothing && return nothing
    return parse(Float64, m.captures[1])
end

"""
Lê um arquivo SDP JSON e retorna (best_mean_fitness, mass_at_best_mean_fitness).
Para cada massa, calcula a média de todos os valores de S (sobre patches, bins
de gordura e tempos). Em seguida, retorna a massa cuja média foi maior.
"""
function best_mean_fitness_and_mass(sdp_file)
    data   = JSON3.read(read(sdp_file, String), Dict)
    S_dict = data["S"]

    best_mean = -Inf
    best_mass = NaN

    for (mass_str, S_nested) in S_dict
        hmax   = length(S_nested)
        f_bins = length(S_nested[1])
        tmax   = length(S_nested[1][1])

        S = Array{Float64}(undef, hmax, f_bins, tmax)
        for h in 1:hmax, y in 1:f_bins, t in 1:tmax
            S[h, y, t] = S_nested[h][y][t]
        end

        fit = mean(S)

        if fit > best_mean
            best_mean = fit
            best_mass = parse(Float64, string(mass_str))
        end
    end

    return best_mean, best_mass
end

"""
Carrega todos os arquivos de um fisio_type e retorna vetores ordenados por mu:
  (mus, max_fitnesses, masses_at_max)
"""
function load_mu_results(fisio_type)
    files = filter(
        f -> occursin(Regex("SDP_$(fisio_type)_mu.+\\.json"), basename(f)),
        readdir(data_dir, join=true)
    )

    if isempty(files)
        @warn "Nenhum arquivo encontrado para '$fisio_type' em $data_dir"
        return Float64[], Float64[], Float64[]
    end

    rows = Tuple{Float64,Float64,Float64}[]
    for fp in sort(files)
        mu = parse_mu(fp, fisio_type)
        mu === nothing && continue
        best_mean, mass = best_mean_fitness_and_mass(fp)
        push!(rows, (mu, best_mean, mass))
        println("  $fisio_type  mu=$mu  mean_fitness=$best_mean  massa=$mass")
    end

    sort!(rows, by = r -> r[1])
    mus    = [r[1] for r in rows]
    fits   = [r[2] for r in rows]
    masses = [r[3] for r in rows]
    return mus, fits, masses
end

# ── Coleta de dados ───────────────────────────────────────────────────────────

println("\n=== Carregando dados ===")
mus_b, fits_b, masses_b = load_mu_results(fisio_type_browser)
mus_g, fits_g, masses_g = load_mu_results(fisio_type_grazer)
mus_m, fits_m, masses_m = load_mu_results(fisio_type_mixed)

# ── Configuração global de estilo ─────────────────────────────────────────────

default(
    fontfamily = "Helvetica",
    guidefont  = font(18),
    tickfont   = font(15),
    legendfont = font(13),
    dpi        = 600
)

browser_color = "#1B9E77"
grazer_color  = "#D95F02"
mixed_color   = "#7570B3"

# ── Plot 1: Max mean fitness vs mu ────────────────────────────────────────────

println("\n=== Gerando gráficos ===")

p1 = plot(mus_g, fits_g,
          xlabel     = "μ",
          ylabel     = "Max mean fitness",
          color      = grazer_color,
          framestyle = :box,
          marker     = :none,
          legend     = :none,
          lw         = 5,
          size       = (900, 700),
          grid       = false)

plot!(p1, mus_b, fits_b, color=browser_color, marker=:none, lw=5)
plot!(p1, mus_m, fits_m, color=mixed_color,   marker=:none, lw=5)

# ── Plot 2: Massa no fitness máximo vs mu ─────────────────────────────────────

yticks_values = [0, 1, 2, 3, 4]
yticks_labels = ["10⁰", "10¹", "10²", "10³", "10⁴"]

p2 = plot(mus_g, log10.(masses_g),
          xlabel     = "μ",
          ylabel     = "Mass at max mean fitness (kg)",
          color      = grazer_color,
          framestyle = :box,
          marker     = :none,
          legend     = :none,
          lw         = 5,
          size       = (900, 700),
          grid       = false,
          yticks     = (yticks_values, yticks_labels),
          ylim       = (0, 4.1))

plot!(p2, mus_b, log10.(masses_b), color=browser_color, marker=:none, lw=5)
plot!(p2, mus_m, log10.(masses_m), color=mixed_color,   marker=:none, lw=5)

# ── Salva figuras ─────────────────────────────────────────────────────────────

results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

savefig(p1, joinpath(results_dir, "max_mean_fitness_vs_mu.svg"))
savefig(p2, joinpath(results_dir, "mass_at_max_mean_fitness_vs_mu.svg"))

println("✓ Gráficos salvos em: $results_dir")