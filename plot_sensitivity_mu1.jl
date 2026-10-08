#!/usr/bin/env julia
"""
plot_sensitivity_mu.jl
----------------------
Lê os arquivos SDP_<fisio_type>_mu<valor>.json gerados pelo SDP.jl, obtém para
cada mu o fitness máximo (F*) e a massa nesse fitness máximo (m*), e plota a
sensibilidade

    ΔF* / Δm*   em função de mu

calculada por diferenças finitas entre valores consecutivos de mu
(o ponto é plotado no mu médio de cada par).

Uso:
  julia plot_sensitivity_mu.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed>
  Exemplo: julia plot_sensitivity_mu.jl browsers grazers mixed
"""

using JSON3
using Statistics
using StatsPlots
using Plots

# ── Opções ────────────────────────────────────────────────────────────────────

# false: Δm em kg (literal ΔF/Δm)
# true : Δ log10(m), coerente com o eixo log do gráfico de massa
const USE_LOG_MASS = false

# compara mu_i com mu_(i+STEP). STEP = 1 usa vizinhos imediatos; aumente
# (2, 3, ...) se m* quase não muda entre vizinhos.
const STEP = 1

# ── Argumentos ────────────────────────────────────────────────────────────────

if length(ARGS) < 3
    println("Uso: julia plot_sensitivity_mu1.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed>")
    println("Exemplo: julia plot_sensitivity_mu1.jl browsers grazers mixed")
    exit(1)
end

fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]

script_dir = @__DIR__
data_dir   = joinpath(script_dir, "data_1")

# ── Funções auxiliares (idênticas ao plot_fitness_mu.jl) ──────────────────────

"""Extrai o valor numérico de mu do nome do arquivo."""
function parse_mu(filename, fisio_type)
    pattern = Regex("SDP_$(fisio_type)_mu(.+)\\.json")
    m = match(pattern, basename(filename))
    m === nothing && return nothing
    return parse(Float64, m.captures[1])
end

"""
Lê um arquivo SDP JSON e retorna (best_mean_fitness, mass_at_best_mean_fitness).
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

"""Retorna (mus, max_fitnesses, masses_at_max) ordenados por mu."""
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
    return [r[1] for r in rows], [r[2] for r in rows], [r[3] for r in rows]
end

"""
Sensibilidade ΔF*/Δm* entre valores consecutivos de mu.
Retorna (mu_medio, sensibilidade). Onde Δm* == 0 (massa ótima não mudou entre
dois mu vizinhos, comum porque as massas estão numa grade discreta), a razão
é indefinida e o ponto vira NaN (a linha fica interrompida ali).
"""
function sensitivity(mus, fits, masses; step = STEP)
    n = length(mus)
    n <= step && return Float64[], Float64[]

    m = USE_LOG_MASS ? log10.(masses) : masses

    idx = 1:(n - step)
    dF  = fits[idx .+ step] .- fits[idx]
    dm  = m[idx .+ step]    .- m[idx]

    sens   = [iszero(d) ? NaN : f / d for (f, d) in zip(dF, dm)]
    mu_mid = (mus[idx] .+ mus[idx .+ step]) ./ 2

    # descarta pontos indefinidos (Δm = 0) para a linha conectar os válidos
    ok = isfinite.(sens)
    println("  pontos válidos: $(count(ok)) de $(length(ok))  " *
            "(massas distintas: $(length(unique(masses))))")
    return mu_mid[ok], sens[ok]
end

# ── Coleta de dados ───────────────────────────────────────────────────────────

println("\n=== Carregando dados ===")
mus_b, fits_b, masses_b = load_mu_results(fisio_type_browser)
mus_g, fits_g, masses_g = load_mu_results(fisio_type_grazer)
mus_m, fits_m, masses_m = load_mu_results(fisio_type_mixed)

mu_sb, sens_b = sensitivity(mus_b, fits_b, masses_b)
mu_sg, sens_g = sensitivity(mus_g, fits_g, masses_g)
mu_sm, sens_m = sensitivity(mus_m, fits_m, masses_m)

# ── Estilo ────────────────────────────────────────────────────────────────────

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

# ── Plot: ΔFitness/ΔMass vs mu ────────────────────────────────────────────────

println("\n=== Gerando gráfico ===")

ylab = USE_LOG_MASS ? "ΔFitness / Δlog₁₀(Mass)" : "ΔFitness / ΔMass (kg⁻¹)"

p = plot(mu_sg, sens_g,
         xlabel     = "μ",
         ylabel     = ylab,
         color      = grazer_color,
         framestyle = :box,
         marker     = :circle,
         markersize = 6,
         legend     = :none,
         lw         = 5,
         size       = (900, 700),
         grid       = false)

plot!(p, mu_sb, sens_b, color=browser_color, marker=:circle, markersize=6, lw=5)
plot!(p, mu_sm, sens_m, color=mixed_color,   marker=:circle, markersize=6, lw=5)
hline!(p, [0], color=:black, ls=:dash, lw=1)

# ── Salva ─────────────────────────────────────────────────────────────────────

results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

savefig(p, joinpath(results_dir, "sensitivity_dfitness_dmass_vs_mu.svg"))

println("✓ Gráfico salvo em: $results_dir")