using JSON3
using Plots
using Statistics

# Função para processar um arquivo JSON e calcular ganhos líquidos
function process_json(filepath)
    data = JSON3.read(read(filepath, String))
    
    results = Dict(
        "masses" => Float64[],
        "habitat_1_net_gains" => Float64[],
        "habitat_2_net_gains" => Float64[]
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
            prob_1 = collect(habitat_data_1[:prob])
            
            # Converter para matriz
            n_bins = length(gains_1)
            prob_matrix_1 = reshape(prob_1, (n_bins, n_bins))
            
            # Ganho líquido ponderado por probabilidade
            weighted_net_gain_1 = 0.0
            total_prob_1 = 0.0
            
            for i in 1:length(gains_1)
                for j in 1:length(costs_1)
                    if prob_matrix_1[i, j] > 0
                        net_gain = gains_1[i] - costs_1[j]
                        weighted_net_gain_1 += net_gain * prob_matrix_1[i, j]
                        total_prob_1 += prob_matrix_1[i, j]
                    end
                end
            end
            
            # Normalizar pela probabilidade total (caso não seja 1.0)
            mean_net_gain_1 = weighted_net_gain_1 / total_prob_1
            
            push!(results["habitat_1_net_gains"], mean_net_gain_1)
            
            println("  Massa $(round(mass, digits=2)) kg | Habitat Fechado: Líquido=$(round(mean_net_gain_1, digits=1)) kcal (prob total=$(round(total_prob_1, digits=3)))")
        else
            push!(results["habitat_1_net_gains"], NaN)
        end
        
        # Habitat 2 (aberto)
        habitat_data_2 = nothing
        if haskey(data[mass_sym], 2) || haskey(data[mass_sym], Symbol("2"))
            habitat_data_2 = get(data[mass_sym], 2, get(data[mass_sym], Symbol("2"), nothing))
        end
        
        if !isnothing(habitat_data_2)
            gains_2 = collect(habitat_data_2[:gains])
            costs_2 = collect(habitat_data_2[:costs])
            prob_2 = collect(habitat_data_2[:prob])
            
            # Converter para matriz
            n_bins = length(gains_2)
            prob_matrix_2 = reshape(prob_2, (n_bins, n_bins))
            
            # Ganho líquido ponderado por probabilidade
            weighted_net_gain_2 = 0.0
            total_prob_2 = 0.0
            
            for i in 1:length(gains_2)
                for j in 1:length(costs_2)
                    if prob_matrix_2[i, j] > 0
                        net_gain = gains_2[i] - costs_2[j]
                        weighted_net_gain_2 += net_gain * prob_matrix_2[i, j]
                        total_prob_2 += prob_matrix_2[i, j]
                    end
                end
            end
            
            # Normalizar pela probabilidade total (caso não seja 1.0)
            mean_net_gain_2 = weighted_net_gain_2 / total_prob_2
            
            push!(results["habitat_2_net_gains"], mean_net_gain_2)
            
            println("  Massa $(round(mass, digits=2)) kg | Habitat Aberto: Líquido=$(round(mean_net_gain_2, digits=1)) kcal (prob total=$(round(total_prob_2, digits=3)))")
        else
            push!(results["habitat_2_net_gains"], NaN)
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
println("CARREGANDO DADOS - GANHOS LÍQUIDOS")
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
    guidefont = font(18),
    tickfont = font(15),
    legendfont = font(13),
    titlefont = font(18),
    dpi = 600
)

# Função para converter número em superscript Unicode
function to_superscript(n)
    superscripts = Dict('0'=>'⁰', '1'=>'¹', '2'=>'²', '3'=>'³', '4'=>'⁴', 
                       '5'=>'⁵', '6'=>'⁶', '7'=>'⁷', '8'=>'⁸', '9'=>'⁹', '-'=>'⁻')
    return join([superscripts[c] for c in string(n)])
end

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
xticks_labels = [i == 0 ? "10⁰" : "10$(to_superscript(i))" for i in xticks_values]

# Função symlog (symmetric log)
# Linear perto de zero, logarítmico para valores grandes
function symlog(x, linthresh=1.0)
    if abs(x) <= linthresh
        return x
    else
        return sign(x) * (linthresh + log10(abs(x) / linthresh))
    end
end

# Função inversa para criar ticks apropriados
function symlog_ticks(y_min, y_max, linthresh=1.0)
    ticks = Float64[]
    labels = String[]
    
    # Determinar ordens de magnitude
    min_order = floor(Int, log10(max(abs(y_min), linthresh)))
    max_order = ceil(Int, log10(max(abs(y_max), linthresh)))
    
    # Adicionar valores negativos
    if y_min < -linthresh
        for order in max_order:-1:1
            val = -10.0^order
            tick_val = symlog(val, linthresh)
            if tick_val >= y_min
                push!(ticks, tick_val)
                push!(labels, "-10$(to_superscript(order))")
            end
        end
    end
    
    # Adicionar região linear (perto de zero)
    if y_min <= 0 && y_max >= 0
        push!(ticks, 0.0)
        push!(labels, "0")
    end
    
    # Adicionar valores positivos
    if y_max > linthresh
        for order in 1:max_order
            val = 10.0^order
            tick_val = symlog(val, linthresh)
            if tick_val <= y_max
                push!(ticks, tick_val)
                push!(labels, "10$(to_superscript(order))")
            end
        end
    end
    
    return (ticks, labels)
end

# Aplicar symlog aos ganhos líquidos
linthresh = 1.0  # Threshold para região linear
for (label, results) in all_results
    results["habitat_1_net_gains_symlog"] = [symlog(g, linthresh) for g in results["habitat_1_net_gains"]]
    results["habitat_2_net_gains_symlog"] = [symlog(g, linthresh) for g in results["habitat_2_net_gains"]]
end

# Encontrar range do eixo Y (em symlog) para cada habitat
y_range_1_symlog = []
y_range_2_symlog = []
for (label, results) in all_results
    append!(y_range_1_symlog, filter(!isnan, results["habitat_1_net_gains_symlog"]))
    append!(y_range_2_symlog, filter(!isnan, results["habitat_2_net_gains_symlog"]))
end

y_min_1_symlog = minimum(y_range_1_symlog)
y_max_1_symlog = maximum(y_range_1_symlog)
y_min_2_symlog = minimum(y_range_2_symlog)
y_max_2_symlog = maximum(y_range_2_symlog)

# Criar ticks para symlog
yticks_1 = symlog_ticks(y_min_1_symlog, y_max_1_symlog, linthresh)
yticks_2 = symlog_ticks(y_min_2_symlog, y_max_2_symlog, linthresh)

# ========== PLOTS ==========

# Plot 1: Ganhos Líquidos - Habitat Fechado
p1 = plot(
    xlabel = "Body mass (kg)",
    ylabel = "Mean net gain (kcal)",
    #title = "Closed Habitat",
    legend = :none,
    size = (900, 700),
    grid = false,
    framestyle = :box,
    xticks = (xticks_values, xticks_labels),
    ylims = (y_min_1_symlog, y_max_1_symlog),
    yticks = yticks_1
    # left_margin = 15Plots.mm,
    # bottom_margin = 10Plots.mm,
    # right_margin = 5Plots.mm,
    # top_margin = 5Plots.mm
)

for (label, results) in all_results
    plot!(p1, masses_log10[label], results["habitat_1_net_gains_symlog"],
          label = label,
          color = colors[label],
          lw = 5)
end

# Linha de break-even (ganho líquido = 0)
hline!(p1, [0], color = :red, linestyle = :dash, linewidth = 2, label = "Break-even")

# Plot 2: Ganhos Líquidos - Habitat Aberto
p2 = plot(
    xlabel = "Body mass (kg)",
    ylabel = "Mean net gain (kcal)",
   #title = "Open Habitat",
    legend = :none,
    size = (900, 700),
    grid = false,
    framestyle = :box,
    xticks = (xticks_values, xticks_labels),
    ylims = (y_min_2_symlog, y_max_2_symlog),
    yticks = yticks_2
    # left_margin = 15Plots.mm,
    # bottom_margin = 10Plots.mm,
    # right_margin = 5Plots.mm,
    # top_margin = 5Plots.mm
)

for (label, results) in all_results
    plot!(p2, masses_log10[label], results["habitat_2_net_gains_symlog"],
          label = label,
          color = colors[label],
          lw = 5)
end

# Linha de break-even (ganho líquido = 0)
hline!(p2, [0], color = :red, linestyle = :dash, linewidth = 2, label = "Break-even")

# ========== SALVAR PLOTS ==========
results_dir = joinpath(current_dir, "results")
isdir(results_dir) || mkpath(results_dir)

# Salvar plots individuais
savefig(p1, joinpath(results_dir, "net_gains_closed_habitat.svg"))
println("\n✓ Plot salvo: net_gains_closed_habitat.svg")

savefig(p2, joinpath(results_dir, "net_gains_open_habitat.svg"))
println("✓ Plot salvo: net_gains_open_habitat.svg")

# Combinar ambos os plots
# plot_combined = plot(p1, p2, 
#                      layout=(1, 2), 
#                      size=(1800, 700),
#                      left_margin = 12Plots.mm,
#                      right_margin = 10Plots.mm,
#                      bottom_margin = 10Plots.mm,
#                      top_margin = 5Plots.mm)

# savefig(plot_combined, joinpath(results_dir, "net_gains_combined.png"))
# println("✓ Plot combinado salvo: net_gains_combined.png")

# Mostrar plot combinado
#display(plot_combined)

println("\n✓ Análise de ganhos líquidos concluída!")
println("="^60)