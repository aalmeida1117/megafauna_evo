#!/usr/bin/env julia
"""
plot_combined.jl
----------------
Gera figura combinada com dois painéis lado a lado:
  Esquerda: Densidade da massa corporal por estratégia alimentar
  Direita:  Massa no fitness máximo vs μ (por estratégia)
Eixo Y compartilhado: Body Mass (kg) em escala log₁₀
"""

using JSON3, Statistics, StatsPlots, Plots, DataFrames, XLSX, KernelDensity

# ── Cores ─────────────────────────────────────────────────────────────────────
browser_color = "#1B9E77"
grazer_color  = "#D95F02"
mixed_color   = "#7570B3"

# ── Estilo global ─────────────────────────────────────────────────────────────
default(
    fontfamily = "Helvetica",
    guidefont  = font(18),
    tickfont   = font(15),
    legendfont = font(13),
    dpi        = 600
)

yticks_values = [0, 1, 2, 3, 4]
yticks_labels = ["10⁰", "10¹", "10²", "10³", "10⁴"]
ylims_shared  = (0, 4.5)

# ══════════════════════════════════════════════════════════════════════════════
# PAINEL ESQUERDO — Densidade da massa corporal
# ══════════════════════════════════════════════════════════════════════════════

df = DataFrame(XLSX.readtable("pnas201907998_s5_movfsn.xlsx", "Sheet 1"; first_row=2))
df.strategy = similar(df.d13C, String)
for i in 1:nrow(df)
    if -25 <= df.d13C[i] < -8
        df.strategy[i] = "browser"
    elseif -8 <= df.d13C[i] < 1
        df.strategy[i] = "mixed"
    elseif df.d13C[i] >= 1
        df.strategy[i] = "grazer"
    else
        df.strategy[i] = "unknown"
    end
end

df.mass_kg  = df.Mass ./ 1000
df.logmass  = log10.(df.mass_kg)

grazers = df[df.strategy .== "grazer",  :logmass]
browsers = df[df.strategy .== "browser", :logmass]
mixed    = df[df.strategy .== "mixed",   :logmass]

# Calcula as densidades manualmente para rotacionar os eixos (x→density, y→logmass)
function kde_xy(data; npoints=300)
    k = kde(data)
    xs = range(minimum(k.x), maximum(k.x), length=npoints)
    ys = pdf(k, collect(xs))
    return collect(xs), ys
end

gx, gy = kde_xy(grazers)
bx, by = kde_xy(browsers)
mx, my = kde_xy(mixed)

# Plot rotacionado: densidade no eixo X, logmassa no eixo Y
p_left = plot(
    gy, gx;                          # densidade=X, logmassa=Y
    xlabel     = "Density",
    ylabel     = "Body Mass (kg)",
    color      = grazer_color,
    lw         = 5,
    legend     = :none,
    framestyle = :box,
    grid       = false,
    yticks     = (yticks_values, yticks_labels),
    ylim       = ylims_shared,
    size       = (700, 650),
    xflip      = true,               # densidade cresce da direita para esquerda
)
plot!(p_left, by, bx; color=browser_color, lw=5)
plot!(p_left, my, mx; color=mixed_color,   lw=5)

# Linhas tracejadas horizontais nas modas
hline!(p_left, [gx[argmax(gy)]]; lw=2, linestyle=:dash, color=grazer_color)
hline!(p_left, [bx[argmax(by)]]; lw=2, linestyle=:dash, color=browser_color)
hline!(p_left, [mx[argmax(my)]]; lw=2, linestyle=:dash, color=mixed_color)

# ══════════════════════════════════════════════════════════════════════════════
# PAINEL DIREITO — Massa no fitness máximo vs μ
# ══════════════════════════════════════════════════════════════════════════════

script_dir = @__DIR__
data_dir   = joinpath(script_dir, "data_1")

if length(ARGS) < 3
    println("Uso: julia plot_combined.jl <fisio_browser> <fisio_grazer> <fisio_mixed>")
    exit(1)
end
fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]

function parse_mu(filename, fisio_type)
    m = match(Regex("SDP_$(fisio_type)_mu(.+)\\.json"), basename(filename))
    m === nothing && return nothing
    return parse(Float64, m.captures[1])
end

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

function load_mu_results(fisio_type)
    files = filter(
        f -> occursin(Regex("SDP_$(fisio_type)_mu.+\\.json"), basename(f)),
        readdir(data_dir, join=true)
    )
    isempty(files) && (@warn "Sem arquivos para '$fisio_type'"; return Float64[], Float64[], Float64[])
    rows = Tuple{Float64,Float64,Float64}[]
    for fp in sort(files)
        mu = parse_mu(fp, fisio_type)
        mu === nothing && continue
        best_mean, mass = best_mean_fitness_and_mass(fp)
        push!(rows, (mu, best_mean, mass))
    end
    sort!(rows, by = r -> r[1])
    return [r[1] for r in rows], [r[2] for r in rows], [r[3] for r in rows]
end

mus_b, _, masses_b = load_mu_results(fisio_type_browser)
mus_g, _, masses_g = load_mu_results(fisio_type_grazer)
mus_m, _, masses_m = load_mu_results(fisio_type_mixed)

p_right = plot(
    mus_g, log10.(masses_g);
    xlabel     = "μ",
    ylabel     = "",           # eixo Y compartilhado — sem rótulo repetido
    color      = grazer_color,
    framestyle = :box,
    marker     = :none,
    legend     = :none,
    lw         = 5,
    grid       = false,
    yticks     = (yticks_values, yticks_labels),
    ylim       = ylims_shared,
    yformatter = _ -> "",      # oculta os tick labels (painel esquerdo já tem)
    size       = (700, 650),
)
plot!(p_right, mus_b, log10.(masses_b); color=browser_color, lw=5)
plot!(p_right, mus_m, log10.(masses_m); color=mixed_color,   lw=5)

# ══════════════════════════════════════════════════════════════════════════════
# COMBINAÇÃO E SALVAMENTO
# ══════════════════════════════════════════════════════════════════════════════

p_combined = plot(p_left, p_right;
    layout = (1, 2),
    size   = (1400, 650),
    left_margin  = 15Plots.mm,
    right_margin = 5Plots.mm,
    bottom_margin = 10Plots.mm,
)

results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)
savefig(p_combined, joinpath(results_dir, "combined_mass_mu.svg"))
println("✓ Figura salva em: $(joinpath(results_dir, "combined_mass_mu.svg"))")