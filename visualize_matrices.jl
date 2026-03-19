using JSON3
using Plots
using Statistics

# ========== PARÂMETROS ==========
if length(ARGS) < 1
    println("Uso: julia visualize_matrices.jl <tipo_fisio> [massa_especifica] [habitat]")
    println("Ex: julia visualize_matrices.jl grazers")
    println("Ex: julia visualize_matrices.jl grazers 10.0 1")
    println("\ntipos disponíveis: grazers, browsers, mixed")
    println("habitat: 1 (fechado) ou 2 (aberto)")
    exit(1)
end

fisio_type = ARGS[1]
massa_especifica = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : nothing
habitat_especifico = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : nothing

# ========== CARREGAR DADOS ==========
current_dir = @__DIR__
json_file = joinpath(current_dir, "data", "$(fisio_type)_dict.json")

if !isfile(json_file)
    error("Arquivo não encontrado: $json_file")
end

println("Carregando: $json_file")
data = JSON3.read(read(json_file, String))

# ========== LISTAR MASSAS DISPONÍVEIS ==========
mass_keys = sort([parse(Float64, string(k)) for k in keys(data)])
println("\n=== MASSAS DISPONÍVEIS ===")
for (i, m) in enumerate(mass_keys)
    println("  $i. $(round(m, digits=2)) kg")
end

# ========== FUNÇÃO PARA VISUALIZAR UMA MATRIZ ==========
function visualize_matrix(mass, habitat, data)
    mass_str = string(mass)
    mass_sym = Symbol(mass_str)
    
    # Acessar dados
    habitat_data = nothing
    if haskey(data[mass_sym], habitat)
        habitat_data = data[mass_sym][habitat]
    elseif haskey(data[mass_sym], Symbol(string(habitat)))
        habitat_data = data[mass_sym][Symbol(string(habitat))]
    else
        println("❌ Dados não encontrados para massa=$mass, habitat=$habitat")
        return
    end
    
    # Extrair dados
    gains = collect(habitat_data[:gains])
    costs = collect(habitat_data[:costs])
    prob = collect(habitat_data[:prob])
    
    # Converter para matriz
    n_bins = length(gains)
    prob_matrix = reshape(prob, (n_bins, n_bins))
    
    # Informações estatísticas
    habitat_name = habitat == 1 ? "Fechado" : "Aberto"
    println("\n" * "="^60)
    println("MASSA: $(round(mass, digits=2)) kg | HABITAT: $habitat_name")
    println("="^60)
    println("Dimensões da matriz: $(size(prob_matrix))")
    println("Número de bins: $n_bins")
    println("\n--- GANHOS ---")
    println("  Min: $(round(minimum(gains), digits=2)) kcal")
    println("  Max: $(round(maximum(gains), digits=2)) kcal")
    println("  Média: $(round(mean(gains), digits=2)) kcal")
    println("\n--- CUSTOS ---")
    println("  Min: $(round(minimum(costs), digits=2)) kcal")
    println("  Max: $(round(maximum(costs), digits=2)) kcal")
    println("  Média: $(round(mean(costs), digits=2)) kcal")
    println("\n--- MATRIZ DE PROBABILIDADES ---")
    println("  Soma total: $(round(sum(prob_matrix), digits=4))")
    println("  Valor máximo: $(round(maximum(prob_matrix), digits=6))")
    println("  Valor mínimo: $(round(minimum(prob_matrix), digits=6))")
    println("  Número de zeros: $(sum(prob_matrix .== 0))")
    
    # Criar heatmap
    p1 = heatmap(
        costs, gains, prob_matrix,
        xlabel="Custos (kcal)",
        ylabel="Ganhos (kcal)",
        title="Matriz de Probabilidades Conjunta\nMassa: $(round(mass, digits=2)) kg | Habitat: $habitat_name",
        color=:viridis,
        size=(800, 700),
        margin=5Plots.mm
    )
    
    # Histograma marginal de ganhos
    gains_marginal = sum(prob_matrix, dims=2) |> vec
    p2 = plot(
        gains, gains_marginal,
        xlabel="Ganhos (kcal)",
        ylabel="Probabilidade Marginal",
        title="Distribuição Marginal de Ganhos",
        linewidth=2,
        color=:green,
        legend=false,
        size=(800, 400)
    )
    
    # Histograma marginal de custos
    costs_marginal = sum(prob_matrix, dims=1) |> vec
    p3 = plot(
        costs, costs_marginal,
        xlabel="Custos (kcal)",
        ylabel="Probabilidade Marginal",
        title="Distribuição Marginal de Custos",
        linewidth=2,
        color=:red,
        legend=false,
        size=(800, 400)
    )
    
    # Scatter plot de ganho vs custo (colorido por probabilidade)
    # Criar pontos apenas onde prob > 0 para visualizar melhor
    idx_nonzero = findall(prob_matrix .> 0)
    gains_points = [gains[i[1]] for i in idx_nonzero]
    costs_points = [costs[i[2]] for i in idx_nonzero]
    prob_points = [prob_matrix[i] for i in idx_nonzero]
    
    p4 = scatter(
        costs_points, gains_points,
        marker_z=prob_points,
        xlabel="Custos (kcal)",
        ylabel="Ganhos (kcal)",
        title="Ganhos vs Custos (colorido por probabilidade)",
        color=:viridis,
        markersize=4,
        alpha=0.6,
        size=(800, 700),
        colorbar_title="Probabilidade",
        legend=false
    )
    
    # Linha de break-even (ganho = custo)
    max_val = max(maximum(gains), maximum(costs))
    plot!(p4, [0, max_val], [0, max_val], 
          linestyle=:dash, color=:white, linewidth=2, label="Break-even")
    
    # Combinar plots
    plot_combined = plot(p1, p2, p3, p4, 
                         layout=(2, 2), 
                         size=(1600, 1400),
                         margin=5Plots.mm)
    
    # Salvar
    output_file = joinpath(current_dir, "matrix_$(fisio_type)_mass$(round(mass, digits=2))_habitat$(habitat).png")
    savefig(plot_combined, output_file)
    println("\n✓ Plot salvo: $output_file")
    
    display(plot_combined)
    
    return prob_matrix, gains, costs
end

# ========== VISUALIZAR ==========
if !isnothing(massa_especifica)
    # Visualizar massa específica
    habitats = !isnothing(habitat_especifico) ? [habitat_especifico] : [1, 2]
    
    for h in habitats
        visualize_matrix(massa_especifica, h, data)
    end
else
    # Mostrar resumo de todas as massas
    println("\n" * "="^60)
    println("RESUMO DE TODAS AS MASSAS - $(uppercase(fisio_type))")
    println("="^60)
    
    for mass in mass_keys
        mass_sym = Symbol(string(mass))
        
        println("\n--- Massa: $(round(mass, digits=2)) kg ---")
        
        for habitat in [1, 2]
            habitat_name = habitat == 1 ? "Fechado" : "Aberto"
            
            # Tentar acessar dados
            habitat_data = nothing
            if haskey(data[mass_sym], habitat)
                habitat_data = data[mass_sym][habitat]
            elseif haskey(data[mass_sym], Symbol(string(habitat)))
                habitat_data = data[mass_sym][Symbol(string(habitat))]
            end
            
            if !isnothing(habitat_data)
                gains = collect(habitat_data[:gains])
                costs = collect(habitat_data[:costs])
                
                println("  $habitat_name:")
                println("    Ganhos: $(round(minimum(gains), digits=1)) - $(round(maximum(gains), digits=1)) kcal (média: $(round(mean(gains), digits=1)))")
                println("    Custos: $(round(minimum(costs), digits=1)) - $(round(maximum(costs), digits=1)) kcal (média: $(round(mean(costs), digits=1)))")
            end
        end
    end
    
    println("\n" * "="^60)
    println("Para visualizar uma matriz específica, use:")
    println("  julia visualize_matrices.jl $(fisio_type) <massa> <habitat>")
    println("Exemplo:")
    println("  julia visualize_matrices.jl $(fisio_type) $(mass_keys[1]) 1")
    println("="^60)
end