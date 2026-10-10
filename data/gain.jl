using JSON3
using Plots
using Statistics

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
    mass_keys = sort([parse(Float64, string(k)) for k in keys(data)])
    
    println("\nProcessando $(length(mass_keys)) massas...")
    
    for mass in mass_keys
        mass_str = string(mass)
        mass_sym = Symbol(mass_str)
        push!(results["masses"], mass)
        
        # Habitat 1 (fechado)
        habitat_data_1 = nothing
        if haskey(data[mass_sym], 1) || haskey(data[mass_sym], Symbol("1"))
            habitat_data_1 = get(data[mass_sym], 1, get(data[mass_sym], Symbol("1"), nothing))
        end
        
        if !isnothing(habitat_data_1)
            gains_1 = collect(habitat_data_1[:gains])
            costs_1 = collect(habitat_data_1[:costs])
            
            # Média dos bins
            mean_gain_1 = mean(gains_1)
            mean_cost_1 = mean(costs_1)
            
            push!(results["habitat_1_gains"], mean_gain_1)
            push!(results["habitat_1_costs"], mean_cost_1)
            
            println("  Massa $(round(mass, digits=2)) kg | Habitat Fechado: Ganho=$(round(mean_gain_1, digits=1)) kcal, Custo=$(round(mean_cost_1, digits=1)) kcal")
        else
            push!(results["habitat_1_gains"], NaN)
            push!(results["habitat_1_costs"], NaN)
        end
        
        # Habitat 2 (aberto)
        habitat_data_2 = nothing
        if haskey(data[mass_sym], 2) || haskey(data[mass_sym], Symbol("2"))
            habitat_data_2 = get(data[mass_sym], 2, get(data[mass_sym], Symbol("2"), nothing))
        end
        
        if !isnothing(habitat_data_2)
            gains_2 = collect(habitat_data_2[:gains])
            costs_2 = collect(habitat_data_2[:costs])
            
            # Média dos bins
            mean_gain_2 = mean(gains_2)
            mean_cost_2 = mean(costs_2)
            
            push!(results["habitat_2_gains"], mean_gain_2)
            push!(results["habitat_2_costs"], mean_cost_2)
            
            println("  Massa $(round(mass, digits=2)) kg | Habitat Aberto: Ganho=$(round(mean_gain_2, digits=1)) kcal, Custo=$(round(mean_cost_2, digits=1)) kcal")
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
    "Browser" => joinpath(current_dir, "data", "browsers_dict.json"),
    "Grazer" => joinpath(current_dir, "data", "grazers_dict.json"),
    "Mixed Feeder" => joinpath(current_dir, "data", "mixed_dict.json")
)

# Processar todos os arquivos
println("="^60)
println("CARREGANDO DADOS")
println("="^60)

all_results = Dict()
for (label, filepath) in json_files
    if isfile(filepath)
        println("\n📁 Processando: $label")
        all_results[label] = process_json(filepath)
    else
        println("\n❌ Arquivo não encontrado: $filepath")
    end
end

# Verificar se há dados para plotar
if isempty(all_results)
    error("Nenhum arquivo JSON foi encontrado!")
end

println("\n" * "="^60)
println("GERANDO PLOTS")
println("="^60)

# ========== CONFIGURAÇÃO GLOBAL DE ESTILO ==========
default(
    fontfamily = "Helvetica",
    guidefont = font(16),
    tickfont = font(13),
    legendfont = font(13),
    titlefont = font(18),
    dpi = 600
)

# Cores do esquema original
browser_color = "#1B9E77"
grazer_color  = "#D95F02"
mixed_color   = "#7570B3"

colors = Dict(
    "Browser" => browser_color,
    "Grazer" => grazer_color,
    "Mixed Feeder" => mixed_color
)

# Converter massas para log10 para o eixo x
masses_log10 = Dict()
for (label, results) in all_results
    masses_log10[label] = log10.(results["masses"])
end

# Definir ticks para o eixo x (potências de 10)
min_exp = minimum([minimum(masses_log10[label]) for label in keys(masses_log10)])
max_exp = maximum([maximum(masses_log10[label]) for label in keys(masses_log10)])
xticks_values = collect(floor(Int, min_exp):ceil(Int, max_exp))
xticks_labels = ["10⁰", "10¹", "10²", "10³", "10⁴"]

# ========== PLOTS ==========

# Plot 1: Ganhos Médios - Habitat Fechado
p1 = plot(
    #xlabel = "Body mass (kg)",
    ylabel = "Mean gain (kcal)",
    title = "Closed Habitat",
    legend = :best,
    size = (800, 650),
    grid = false,
    yscale = :log10,
    xticks = (xticks_values, xticks_labels)
)

for (label, results) in all_results
    plot!(p1, masses_log10[label], results["habitat_1_gains"],
          label = label,
          color = colors[label],
          lw = 2.5)
end

# Plot 2: Ganhos Médios - Habitat Aberto
p2 = plot(
    #xlabel = "Body mass (kg)",
    #ylabel = "Mean gain (kcal)",
    title = "Open Habitat",
    legend = :best,
    size = (800, 650),
    grid = false,
    yscale = :log10,
    xticks = (xticks_values, xticks_labels)
)

for (label, results) in all_results
    plot!(p2, masses_log10[label], results["habitat_2_gains"],
          label = label,
          color = colors[label],
          lw = 2.5)
end

# Plot 3: Custos Médios - Habitat Fechado
p3 = plot(
    xlabel = "Body mass (kg)",
    ylabel = "Mean cost (kcal)",
    #title = "Mean Costs - Closed Habitat",
    legend = :best,
    size = (800, 650),
    grid = false,
    yscale = :log10,
    xticks = (xticks_values, xticks_labels)
)

for (label, results) in all_results
    plot!(p3, masses_log10[label], results["habitat_1_costs"],
          label = label,
          color = colors[label],
          lw = 2.5)
end

# Plot 4: Custos Médios - Habitat Aberto
p4 = plot(
    xlabel = "Body mass (kg)",
    #ylabel = "Mean cost (kcal)",
    #title = "Mean Costs - Open Habitat",
    legend = :best,
    size = (800, 650),
    grid = false,
    yscale = :log10,
    xticks = (xticks_values, xticks_labels)
)

for (label, results) in all_results
    plot!(p4, masses_log10[label], results["habitat_2_costs"],
          label = label,
          color = colors[label],
          lw = 2.5)
end

# ========== SALVAR PLOTS ==========
results_dir = joinpath(current_dir, "results")
isdir(results_dir) || mkpath(results_dir)

# Combinar todos os plots em uma única imagem
plot_combined = plot(p1, p2, p3, p4, 
                     layout=(2, 2), 
                     size=(1800, 1600),
                     left_margin = 15Plots.mm,
                     right_margin = 10Plots.mm,
                     top_margin = 10Plots.mm,
                     bottom_margin = 15Plots.mm)

output_file = joinpath(results_dir, "gains_costs_combined.png")
savefig(plot_combined, output_file)
println("\n✓ Plot combinado salvo: $output_file")

# Mostrar plot combinado
display(plot_combined)

println("\n✓ Análise concluída!")
println("="^60)