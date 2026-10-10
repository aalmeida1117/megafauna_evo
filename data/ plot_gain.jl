using JSON3
using Plots
using Statistics

# Função para calcular média ponderada pela probabilidade
function weighted_mean(values, probabilities)
    return sum(values .* probabilities)
end

# Função para processar um arquivo JSON
function process_json(filepath)
    data = JSON3.read(read(filepath, String))
    
    results = Dict(
        "masses" => Float64[],
        "habitat_1_gains" => Float64[],
        "habitat_1_costs" => Float64[],
        "habitat_2_gains" => Float64[],
        "habitat_2_costs" => Float64[]
    )
    
    # Iterar pelas massas (ordenar para manter consistência)
    mass_keys = sort([parse(Float64, k) for k in keys(data)])
    
    for mass in mass_keys
        mass_str = string(mass)
        push!(results["masses"], mass)
        
        # Habitat 1 (fechado)
        if haskey(data[mass_str], 1)
            gains_1 = data[mass_str][1]["gains"]
            costs_1 = data[mass_str][1]["costs"]
            prob_1 = data[mass_str][1]["prob"]
            
            # Calcular média marginal para ganhos (soma sobre custos)
            gains_marginal = sum(prob_1, dims=2) |> vec
            mean_gain_1 = weighted_mean(gains_1, gains_marginal)
            
            # Calcular média marginal para custos (soma sobre ganhos)
            costs_marginal = sum(prob_1, dims=1) |> vec
            mean_cost_1 = weighted_mean(costs_1, costs_marginal)
            
            push!(results["habitat_1_gains"], mean_gain_1)
            push!(results["habitat_1_costs"], mean_cost_1)
        else
            push!(results["habitat_1_gains"], NaN)
            push!(results["habitat_1_costs"], NaN)
        end
        
        # Habitat 2 (aberto)
        if haskey(data[mass_str], 2)
            gains_2 = data[mass_str][2]["gains"]
            costs_2 = data[mass_str][2]["costs"]
            prob_2 = data[mass_str][2]["prob"]
            
            gains_marginal = sum(prob_2, dims=2) |> vec
            mean_gain_2 = weighted_mean(gains_2, gains_marginal)
            
            costs_marginal = sum(prob_2, dims=1) |> vec
            mean_cost_2 = weighted_mean(costs_2, costs_marginal)
            
            push!(results["habitat_2_gains"], mean_gain_2)
            push!(results["habitat_2_costs"], mean_cost_2)
        else
            push!(results["habitat_2_gains"], NaN)
            push!(results["habitat_2_costs"], NaN)
        end
    end
    
    return results
end

# Obter diretório atual
current_dir = @__DIR__

# Caminhos dos arquivos JSON
json_files = Dict(
    "Grazers" => joinpath(current_dir, "data", "grazers_dict.json"),
    "Browsers" => joinpath(current_dir, "data", "browsers_dict.json"),
    "Mixed Feeders" => joinpath(current_dir, "data", "mixed_dict.json")
)

# Processar todos os arquivos
all_results = Dict()
for (label, filepath) in json_files
    if isfile(filepath)
        println("Processando: $filepath")
        all_results[label] = process_json(filepath)
    else
        println("Arquivo não encontrado: $filepath")
    end
end

# Verificar se há dados para plotar
if isempty(all_results)
    error("Nenhum arquivo JSON foi encontrado!")
end

# ========== PLOTS ==========

# Configuração de estilo
colors = Dict(
    "Grazers" => :green,
    "Browsers" => :brown,
    "Mixed Feeders" => :purple
)

markers = Dict(
    "Grazers" => :circle,
    "Browsers" => :square,
    "Mixed Feeders" => :diamond
)

# Plot 1: Ganhos Médios - Habitat Fechado (log-log)
p1 = plot(
    xlabel="Massa (kg)",
    ylabel="Ganho Médio (kcal)",
    title="Ganhos Médios - Habitat Fechado",
    xscale=:log10,
    yscale=:log10,
    legend=:bottomright,
    size=(800, 600),
    grid=true,
    minorgrid=true
)

for (label, results) in all_results
    plot!(p1, results["masses"], results["habitat_1_gains"],
          label=label,
          color=colors[label],
          marker=markers[label],
          markersize=6,
          linewidth=2,
          alpha=0.8)
end

# Plot 2: Ganhos Médios - Habitat Aberto (log-log)
p2 = plot(
    xlabel="Massa (kg)",
    ylabel="Ganho Médio (kcal)",
    title="Ganhos Médios - Habitat Aberto",
    xscale=:log10,
    yscale=:log10,
    legend=:bottomright,
    size=(800, 600),
    grid=true,
    minorgrid=true
)

for (label, results) in all_results
    plot!(p2, results["masses"], results["habitat_2_gains"],
          label=label,
          color=colors[label],
          marker=markers[label],
          markersize=6,
          linewidth=2,
          alpha=0.8)
end

# Plot 3: Custos Médios - Habitat Fechado (log-log)
p3 = plot(
    xlabel="Massa (kg)",
    ylabel="Custo Médio (kcal)",
    title="Custos Médios - Habitat Fechado",
    xscale=:log10,
    yscale=:log10,
    legend=:bottomright,
    size=(800, 600),
    grid=true,
    minorgrid=true
)

for (label, results) in all_results
    plot!(p3, results["masses"], results["habitat_1_costs"],
          label=label,
          color=colors[label],
          marker=markers[label],
          markersize=6,
          linewidth=2,
          alpha=0.8)
end

# Plot 4: Custos Médios - Habitat Aberto (log-log)
p4 = plot(
    xlabel="Massa (kg)",
    ylabel="Custo Médio (kcal)",
    title="Custos Médios - Habitat Aberto",
    xscale=:log10,
    yscale=:log10,
    legend=:bottomright,
    size=(800, 600),
    grid=true,
    minorgrid=true
)

for (label, results) in all_results
    plot!(p4, results["masses"], results["habitat_2_costs"],
          label=label,
          color=colors[label],
          marker=markers[label],
          markersize=6,
          linewidth=2,
          alpha=0.8)
end

# Plot 5: Comparação Ganho Líquido (Ganho - Custo) - Habitat Fechado
p5 = plot(
    xlabel="Massa (kg)",
    ylabel="Ganho Líquido (kcal)",
    title="Ganho Líquido - Habitat Fechado",
    xscale=:log10,
    legend=:bottomright,
    size=(800, 600),
    grid=true,
    minorgrid=true
)

for (label, results) in all_results
    net_gain = results["habitat_1_gains"] .- results["habitat_1_costs"]
    plot!(p5, results["masses"], net_gain,
          label=label,
          color=colors[label],
          marker=markers[label],
          markersize=6,
          linewidth=2,
          alpha=0.8)
end
hline!(p5, [0], color=:red, linestyle=:dash, label="Break-even", linewidth=1.5)

# Plot 6: Comparação Ganho Líquido - Habitat Aberto
p6 = plot(
    xlabel="Massa (kg)",
    ylabel="Ganho Líquido (kcal)",
    title="Ganho Líquido - Habitat Aberto",
    xscale=:log10,
    legend=:bottomright,
    size=(800, 600),
    grid=true,
    minorgrid=true
)

for (label, results) in all_results
    net_gain = results["habitat_2_gains"] .- results["habitat_2_costs"]
    plot!(p6, results["masses"], net_gain,
          label=label,
          color=colors[label],
          marker=markers[label],
          markersize=6,
          linewidth=2,
          alpha=0.8)
end
hline!(p6, [0], color=:red, linestyle=:dash, label="Break-even", linewidth=1.5)

# Combinar todos os plots
plot_combined = plot(p1, p2, p3, p4, p5, p6, 
                     layout=(3, 2), 
                     size=(1600, 1800),
                     margin=5Plots.mm)

# Salvar
output_file = joinpath(current_dir, "plots_ganhos_custos_completo.png")
savefig(plot_combined, output_file)
println("\nPlot salvo: $output_file")

# Mostrar plots individuais também
display(p1)
display(p2)
display(p5)
display(p6)

println("\n✓ Análise concluída!")