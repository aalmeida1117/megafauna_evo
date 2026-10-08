using JSON3
using Plots

# ============================================================
# Configuração global de estilo (fonte, tamanhos, etc.)
# ============================================================
default(
    fontfamily = "Helvetica",
    guidefont  = font(18),
    tickfont   = font(15),
    legendfont = font(13),
    dpi        = 600
)

# ============================================================
# Paleta de cores (ColorBrewer Dark2)
# ============================================================
browser_color = "#1B9E77"
grazer_color  = "#D95F02"
mixed_color   = "#7570B3"

# ========== CONFIGURAÇÃO ==========
# Cada entrada: (fisio_type, massa EXATAMENTE como aparece nas chaves do JSON, rótulo no gráfico, cor)
# Atenção: as chaves de massa no JSON costumam ter casa decimal (ex: "100.0", "1.0"),
# então massa 10 provavelmente é a chave "10.0" e não "10" — confira no seu arquivo.
series_config = [
    ("grazers", "6.31", "Grazers (mass = 6.31 kg)", grazer_color),
    ("mixed",   "10.0", "Mixed feeders (mass = 10.0 kg)", mixed_color),
    ("browsers", "25.12", "Browsers (mass = 25.12 kg)", browser_color),
]

data_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "data_1")

function mu_vs_fitness(fisio_type::String, target_mass::String, data_dir::String)
    # Nome esperado: SDP_<fisio_type>_mu<valor>.json  (aceita ponto ou underscore no valor de mu)
    files = filter(f -> startswith(f, "SDP_$(fisio_type)_mu") && endswith(f, ".json"), readdir(data_dir))

    if isempty(files)
        println("Nenhum arquivo SDP_$(fisio_type)_mu*.json encontrado em '$data_dir'.")
        return Tuple{Float64,Float64}[]
    end

    mu_fitness = Tuple{Float64,Float64}[]

    for fname in files
        m = match(r"mu(\d+[._]\d+)\.json$", fname)
        if m === nothing
            println("Aviso: não consegui extrair mu de '$fname', pulando.")
            continue
        end
        mu_val = parse(Float64, replace(m.captures[1], "_" => "."))

        filepath = joinpath(data_dir, fname)
        data = JSON3.read(read(filepath, String))
        S = data["S"]

        if !haskey(S, target_mass)
            println("Aviso: massa '$target_mass' não encontrada em '$fname' ($fisio_type), pulando.")
            continue
        end

        # S[target_mass] tem formato [h][f][t] -> média sobre todos os elementos
        S_mass = S[target_mass]
        total = 0.0
        count = 0
        for h in S_mass
            for f in h
                for t in f
                    total += t
                    count += 1
                end
            end
        end

        push!(mu_fitness, (mu_val, total / count))
    end

    sort!(mu_fitness, by = x -> x[1])
    return mu_fitness
end

# ========== CALCULA E PLOTA CADA SÉRIE ==========
p = plot(
    xlabel     = "μ",
    ylabel     = "Mean Fitness",
    #title      = "Fitness médio vs μ",
    lw         = 5, 
    legend     = :bottomright,
    grid       = false,
    framestyle = :box,
    size       = (800, 650),
)

for (fisio_type, target_mass, label, color) in series_config
    mu_fitness = mu_vs_fitness(fisio_type, target_mass, data_dir)

    if isempty(mu_fitness)
        println("Aviso: nenhuma curva gerada para '$label' — verifique os arquivos/chave de massa.")
        continue
    end

    println("\nμ vs Fitness médio — $label:")
    for (mu_val, f) in mu_fitness
        println("  μ = $mu_val  ->  fitness médio = $f")
    end

    mus     = [x[1] for x in mu_fitness]
    fitness = [x[2] for x in mu_fitness]

    plot!(p, mus, fitness, label = label, lw = 5, color = color)
end

results_dir = joinpath(@__DIR__, "results")
isdir(results_dir) || mkpath(results_dir)
output_file = joinpath(results_dir, "fitness_vs_mu_combined.svg")
savefig(p, output_file)

println("\n✓ Gráfico salvo em: $output_file")