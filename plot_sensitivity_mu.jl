#!/usr/bin/env julia
"""
plot_sensitivity_mu.jl
----------------------
Lê os arquivos SDP_<fisio_type>_mu<valor>.json gerados pelo SDP.jl e plota:
  Δ mean fitness / Δ massa vs mu  (uma curva por tipo: grazer, browser, mixed feeder)

Para cada mu:
  1. Calcula o mean fitness de cada massa (média sobre patches, f_bins e tmax)
  2. Identifica a massa com maior mean fitness (índice i*)
  3. Calcula a derivada discreta nesse ponto usando as massas adjacentes:
       - Se i* é interior: (fit[i*+1] - fit[i*-1]) / (mass[i*+1] - mass[i*-1])  [central]
       - Se i* é o primeiro: (fit[i*+1] - fit[i*]) / (mass[i*+1] - mass[i*])    [forward]
       - Se i* é o último:   (fit[i*] - fit[i*-1]) / (mass[i*] - mass[i*-1])    [backward]

Uso:
  julia plot_sensitivity_mu.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed>
  Exemplo: julia plot_sensitivity_mu.jl browsers grazers mixed
"""

using JSON3
using Statistics
using StatsPlots
using Plots

# ── Argumentos ────────────────────────────────────────────────────────────────

if length(ARGS) < 3
    println("Uso: julia plot_sensitivity_mu.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed>")
    println("Exemplo: julia plot_sensitivity_mu.jl browsers grazers mixed")
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
Lê um arquivo SDP JSON e retorna, para cada massa, o mean fitness.
Retorna um vetor de (massa, mean_fitness) ordenado por massa.
"""
function load_mean_fitness_per_mass(sdp_file)
    data   = JSON3.read(read(sdp_file, String), Dict)
    S_dict = data["S"]

    pairs = Tuple{Float64, Float64}[]

    for (mass_str, S_nested) in S_dict
        hmax   = length(S_nested)
        f_bins = length(S_nested[1])
        tmax   = length(S_nested[1][1])

        S = Array{Float64}(undef, hmax, f_bins, tmax)
        for h in 1:hmax, y in 1:f_bins, t in 1:tmax
            S[h, y, t] = S_nested[h][y][t]
        end

        push!(pairs, (parse(Float64, string(mass_str)), mean(S)))
    end

    sort!(pairs, by = p -> p[1])
    return pairs
end

"""
Dado um vetor de (massa, mean_fitness) ordenado por massa,
retorna (max_mean_fitness, massa_no_maximo).
"""
function max_mean_and_mass(pairs)
    i_star = argmax([p[2] for p in pairs])
    return pairs[i_star][2], pairs[i_star][1]
end

"""
Carrega todos os arquivos de um fisio_type e retorna:
  - mus_mid: pontos médios entre mus adjacentes  (comprimento n-1)
  - sens:    Δ max_mean_fitness / Δ mass*  entre mus adjacentes

Para cada par de mus consecutivos (mu_i, mu_{i+1}):
  sens[i] = (max_fit[i+1] - max_fit[i]) / (mass*[i+1] - mass*[i])
plotado no ponto médio (mu_i + mu_{i+1}) / 2.
"""
function load_sensitivity_results(fisio_type)
    files = filter(
        f -> occursin(Regex("SDP_$(fisio_type)_mu.+\\.json"), basename(f)),
        readdir(data_dir, join=true)
    )

    if isempty(files)
        @warn "Nenhum arquivo encontrado para '$fisio_type' em $data_dir"
        return Float64[], Float64[]
    end

    # coleta (mu, max_mean_fitness, mass*) para cada arquivo
    rows = Tuple{Float64, Float64, Float64}[]
    for fp in sort(files)
        mu = parse_mu(fp, fisio_type)
        mu === nothing && continue
        pairs    = load_mean_fitness_per_mass(fp)
        fit, m   = max_mean_and_mass(pairs)
        push!(rows, (mu, fit, m))
        println("  $fisio_type  mu=$mu  max_mean_fit=$fit  massa*=$m kg")
    end

    sort!(rows, by = r -> r[1])
    n = length(rows)

    if n < 2
        @warn "Menos de 2 valores de mu para '$fisio_type' — não é possível calcular diferenças."
        return Float64[], Float64[]
    end

    mus_mid = Float64[]
    sens    = Float64[]
    for i in 1:(n-1)
        mu1, fit1, m1 = rows[i]
        mu2, fit2, m2 = rows[i+1]
        dm = m2 - m1
        if dm == 0.0
            @warn "Δmassa = 0 entre mu=$mu1 e mu=$mu2 para '$fisio_type' — pulando."
            continue
        end
        push!(mus_mid, (mu1 + mu2) / 2)
        push!(sens,    (fit2 - fit1) / dm)
    end

    return mus_mid, sens
end

# ── Coleta de dados ───────────────────────────────────────────────────────────

println("\n=== Carregando dados ===")
mus_b, sens_b = load_sensitivity_results(fisio_type_browser)
mus_g, sens_g = load_sensitivity_results(fisio_type_grazer)
mus_m, sens_m = load_sensitivity_results(fisio_type_mixed)

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

# ── Plot ──────────────────────────────────────────────────────────────────────

println("\n=== Gerando gráfico ===")

p = plot(mus_g, sens_g,
         xlabel     = "μ",
         ylabel     = "Δ max mean fitness / Δ mass (kg⁻¹)",
         color      = grazer_color,
         framestyle = :box,
         marker     = :none,
         legend     = :none,
         lw         = 5,
         size       = (900, 700),
         grid       = false)

plot!(p, mus_b, sens_b, color=browser_color, marker=:none, lw=5)
plot!(p, mus_m, sens_m, color=mixed_color,   marker=:none, lw=5)

# Linha em y=0 para referência
hline!(p, [0.0], lw=2, linestyle=:dash, color=:black)

# ── Salva figura ──────────────────────────────────────────────────────────────

results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

savefig(p, joinpath(results_dir, "sensitivity_vs_mu.svg"))
println("✓ Gráfico salvo em: $(joinpath(results_dir, "sensitivity_vs_mu.svg"))")