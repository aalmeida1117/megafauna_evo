#!/usr/bin/env julia
"""
plot_mass_at_peak_sensitivity_mu.jl
-----------------------------------
Pergunta: a massa em que fica o pico de ΔFitness/ΔMass muda com mu?

Lê os arquivos SDP_<fisio_type>_mu<valor>.json gerados pelo SDP.jl e, para cada mu:
  1. calcula o fitness médio F(m) de todas as massas (média de S);
  2. calcula ΔF/Δm entre massas consecutivas;
  3. encontra a massa em que ΔF/Δm é máximo (pico da sensibilidade).
Plota essa massa em função de mu (uma curva por tipo: grazer, browser, mixed).

A massa associada a cada ΔF/Δm é a média geométrica das duas massas
consecutivas usadas no cálculo.

Uso:
  julia plot_mass_at_peak_sensitivity_mu.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed>
  Exemplo: julia plot_mass_at_peak_sensitivity_mu.jl browsers grazers mixed
"""

using JSON3
using Statistics
using StatsPlots
using Plots

# ── Opções ────────────────────────────────────────────────────────────────────

# false: Δm em kg (literal ΔF/Δm)
# true : Δ log10(m)
# ATENÇÃO: a posição do pico depende dessa escolha (ΔF/Δm e ΔF/Δlog m não têm
# o máximo na mesma massa). Vale testar as duas e ver se a conclusão se mantém.
const USE_LOG_MASS = false

# ── Argumentos ────────────────────────────────────────────────────────────────

if length(ARGS) < 3
    println("Uso: julia plot_mass_at_peak_sensitivity_mu.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed>")
    println("Exemplo: julia plot_mass_at_peak_sensitivity_mu.jl browsers grazers mixed")
    exit(1)
end

fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]

script_dir = @__DIR__
data_dir   = joinpath(script_dir, "data_1/")

# ── Funções auxiliares ────────────────────────────────────────────────────────

"""Extrai o valor numérico de mu do nome do arquivo."""
function parse_mu(filename, fisio_type)
    pattern = Regex("SDP_$(fisio_type)_mu(.+)\\.json")
    m = match(pattern, basename(filename))
    m === nothing && return nothing
    return parse(Float64, m.captures[1])
end

"""Lê um arquivo SDP JSON e retorna Dict massa => fitness médio (média de S)."""
function mean_fitness_by_mass(sdp_file)
    data   = JSON3.read(read(sdp_file, String), Dict)
    S_dict = data["S"]

    out = Dict{Float64,Float64}()

    for (mass_str, S_nested) in S_dict
        hmax   = length(S_nested)
        f_bins = length(S_nested[1])
        tmax   = length(S_nested[1][1])

        S = Array{Float64}(undef, hmax, f_bins, tmax)
        for h in 1:hmax, y in 1:f_bins, t in 1:tmax
            S[h, y, t] = S_nested[h][y][t]
        end

        out[parse(Float64, string(mass_str))] = mean(S)
    end

    return out
end

"""Carrega todos os arquivos de um fisio_type. Retorna (mus, Vector{Dict}) ordenado por mu."""
function load_mu_results(fisio_type)
    files = filter(
        f -> occursin(Regex("SDP_$(fisio_type)_mu.+\\.json"), basename(f)),
        readdir(data_dir, join=true)
    )

    if isempty(files)
        @warn "Nenhum arquivo encontrado para '$fisio_type' em $data_dir"
        return Float64[], Dict{Float64,Float64}[]
    end

    rows = Tuple{Float64,Dict{Float64,Float64}}[]
    for fp in sort(files)
        mu = parse_mu(fp, fisio_type)
        mu === nothing && continue
        d = mean_fitness_by_mass(fp)
        push!(rows, (mu, d))
        println("  $fisio_type  mu=$mu  ($(length(d)) massas)")
    end

    sort!(rows, by = r -> r[1])
    return [r[1] for r in rows], [r[2] for r in rows]
end

"""
Para cada mu, retorna a massa (média geométrica do par de massas consecutivas)
em que ΔF/Δm é máximo, e o valor desse pico.
"""
function mass_at_peak_sensitivity(mus, dicts)
    peak_mass = fill(NaN, length(mus))
    peak_val  = fill(NaN, length(mus))

    for (j, d) in enumerate(dicts)
        masses = sort(collect(keys(d)))
        length(masses) < 2 && continue
        F  = [d[m] for m in masses]
        x  = USE_LOG_MASS ? log10.(masses) : masses

        sens  = diff(F) ./ diff(x)
        m_mid = sqrt.(masses[1:end-1] .* masses[2:end])

        i = argmax(sens)
        peak_mass[j] = m_mid[i]
        peak_val[j]  = sens[i]
    end

    println("  massas de pico distintas: $(length(unique(filter(!isnan, peak_mass)))) " *
            "(em $(length(mus)) valores de mu)")
    return peak_mass, peak_val
end

# ── Coleta de dados ───────────────────────────────────────────────────────────

println("\n=== Carregando dados ===")
mus_b, dicts_b = load_mu_results(fisio_type_browser)
mus_g, dicts_g = load_mu_results(fisio_type_grazer)
mus_m, dicts_m = load_mu_results(fisio_type_mixed)

println("\n=== Calculando massa no pico de ΔF/Δm ===")
peak_b, _ = isempty(mus_b) ? (Float64[], Float64[]) : mass_at_peak_sensitivity(mus_b, dicts_b)
peak_g, _ = isempty(mus_g) ? (Float64[], Float64[]) : mass_at_peak_sensitivity(mus_g, dicts_g)
peak_m, _ = isempty(mus_m) ? (Float64[], Float64[]) : mass_at_peak_sensitivity(mus_m, dicts_m)

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

# ── Plot: massa no pico de ΔF/Δm vs mu ────────────────────────────────────────

println("\n=== Gerando gráfico ===")

yticks_values = [0, 1, 2, 3, 4]
yticks_labels = ["10⁰", "10¹", "10²", "10³", "10⁴"]

p = plot(xlabel     = "μ",
         ylabel     = "Mass at max ΔFitness/ΔMass (kg)",
         framestyle = :box,
         legend     = :none,
         size       = (900, 700),
         grid       = false,
         yticks     = (yticks_values, yticks_labels),
         ylim       = (0, 4.1))

isempty(mus_g) || plot!(p, mus_g, log10.(peak_g), color=grazer_color,  marker=:none, lw=5)
isempty(mus_b) || plot!(p, mus_b, log10.(peak_b), color=browser_color, marker=:none, lw=5)
isempty(mus_m) || plot!(p, mus_m, log10.(peak_m), color=mixed_color,   marker=:none, lw=5)

# ── Salva ─────────────────────────────────────────────────────────────────────

results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

savefig(p, joinpath(results_dir, "mass_at_max_dfitness_dmass_vs_mu.svg"))

println("✓ Gráfico salvo em: $results_dir")