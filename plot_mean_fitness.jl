using JSON3
using Plots

const tmax = 30
const f_bins = 50

function open_files(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    script_dir = @__DIR__
    data_dir = joinpath(script_dir, "data")
    output_file = joinpath(data_dir, "SDP_$(fisio_type)_$(run_date).json")

    if !isfile(output_file)
        println("Erro: Arquivo '$output_file' não encontrado.")
        return Dict{String, Any}()
    end

    println("Carregando dados de: $output_file")
    data = JSON3.read(read(output_file, String))
    Smass = data["S"]

    return Smass
end 

function mean_fitness(Smass, exp_start::Float64, exp_step::Float64, exp_end::Float64)
    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents
    fitness_values = Dict{Float64, Float64}()
    
    for mass in masses
        mass_str = string(round(mass; digits=2))
        if haskey(Smass, mass_str)
            hmax = length(Smass[mass_str])
            total_fitness = 0.0
            count = 0
            
            for h in 1:hmax
                for i in 1:f_bins
                    for t in 1:tmax
                        total_fitness += Smass[mass_str][h][i][t]
                        count += 1
                    end
                end
            end
            
            fitness_values[mass] = count > 0 ? total_fitness / count : 0.0
        end
    end
     
    return fitness_values
end

function get_fitness_arrays(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    Smass = open_files(fisio_type, exp_start, exp_step, exp_end, run_date)
    
    if isempty(Smass)
        return (Float64[], Float64[])
    end
    
    fitness_dict = mean_fitness(Smass, exp_start, exp_step, exp_end)
    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents
    
    fitness_array = [fitness_dict[m] for m in masses]
    
    return (exponents, fitness_array)
end

function print_fitness_table(exponents, fitness_browsers, fitness_grazers, fitness_mixed)
    println("\n" * "="^80)
    println("VALORES DE FITNESS MÉDIO POR MASSA")
    println("="^80)
    println(rpad("Massa (kg)", 15), " | ", 
            rpad("Browser", 15), " | ", 
            rpad("Grazer", 15), " | ", 
            rpad("Mixed Feeder", 15))
    println("-"^80)
    
    for i in 1:length(exponents)
        mass = 10^exponents[i]
        println(rpad(string(mass), 15), " | ", 
                rpad(string(round(fitness_browsers[i], digits=6)), 15), " | ",
                rpad(string(round(fitness_grazers[i], digits=6)), 15), " | ",
                rpad(string(round(fitness_mixed[i], digits=6)), 15))
    end
    println("="^80 * "\n")
end

# ============================================================
# MAIN
# ============================================================
if length(ARGS) < 7
    println("Uso: julia plot_mean_fitness.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed> <exp_start> <exp_end> <exp_step> <run_date>")
    println("Exemplo: julia plot_mean_fitness.jl browsers grazers mixed 0 4 0.1 20261022")
    exit(1)
end

fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]
exp_start = parse(Float64, ARGS[4])
exp_end   = parse(Float64, ARGS[5])
exp_step  = parse(Float64, ARGS[6])
run_date  = ARGS[7]

println("\n=== Carregando fitness médio ===")

# Carrega os dados dos três morfotipos
exponents_b, fitness_browsers = get_fitness_arrays(fisio_type_browser, exp_start, exp_step, exp_end, run_date)
exponents_g, fitness_grazers  = get_fitness_arrays(fisio_type_grazer, exp_start, exp_step, exp_end, run_date)
exponents_m, fitness_mixed    = get_fitness_arrays(fisio_type_mixed, exp_start, exp_step, exp_end, run_date)

# Imprime tabela com os valores
print_fitness_table(exponents_b, fitness_browsers, fitness_grazers, fitness_mixed)

println("\n=== Gerando gráfico de fitness médio ===")

# ============================================================
# Configuração global de estilo (fonte, tamanhos, etc.)
# ============================================================
default(
    fontfamily = "Helvetica",
    guidefont = font(18),
    tickfont = font(15),
    legendfont = font(13),
    dpi = 600
)

# ============================================================
# Plotagem aprimorada
# ============================================================
browser_color = "#1B9E77"
grazer_color  = "#D95F02"
mixed_color   = "#7570B3"

# Define os ticks desejados (somente potências inteiras de 10)
xticks_values = [0, 1, 2, 3, 4]  # eixo está em log10(massa)
xticks_labels = ["10⁰", "10¹", "10²", "10³", "10⁴"]

p = plot(
    exponents_b, fitness_browsers,
    xlabel = "Body mass (kg)",
    ylabel = "Mean fitness",
    lw = 5,
    marker = :none,
    color = browser_color,
    ylim = (0, 1),
    legend = :none,
    size = (800, 650),
    grid = false,
    xlim = (-0.2, 4.2),
    framestyle = :box,
    xticks = (xticks_values, xticks_labels)   # ← só 10⁰–10⁴
)

plot!(p, exponents_g, fitness_grazers,
   #label = "Grazer",
    lw = 5,
    marker = :none,
    color = grazer_color
)

plot!(p, exponents_m, fitness_mixed,
    #label = "Mixed Feeder",
    lw = 5,
    marker = :none,
    color = mixed_color
)

vline!(p, [log10(3.5)],
    lw =2,
    linestyle = :dash,
    color = :red
)

vline!(p, [log10(6.5)],
    lw =2,
    linestyle = :dash,
    color = :red
)

# ============================================================
# Salvando o gráfico
# ============================================================
script_dir = @__DIR__
results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

output_file = joinpath(results_dir, "mean_fitness_$(run_date).svg")
savefig(p, output_file)

println("✓ Gráfico de fitness médio salvo: $(output_file)")
println("\n✅ Análise concluída com sucesso!")