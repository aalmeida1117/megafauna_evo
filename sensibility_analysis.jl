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
        return (Dict{String, Any}(), Dict{String, Any}())
    end

    println("Carregando dados de: $output_file")
    data = JSON3.read(read(output_file, String))
    Tmass = data["T"]
    Smass = data["S"]

    return (Tmass, Smass)
end 

function prop_open_habitat(Tmass, exp_start::Float64, exp_step::Float64, exp_end::Float64)
    proportions = Dict{Float64, Float64}()

    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents
    
    for mass in masses
        mass_str = string(round(mass; digits=2))
        if haskey(Tmass, mass_str)
            hmax_actual = length(Tmass[mass_str])
            T_matrix_for_mass = Array{Float64}(undef, hmax_actual, f_bins, tmax)

            for h_idx in 1:hmax_actual
                for f_idx in 1:f_bins
                    for t_idx in 1:tmax
                        T_matrix_for_mass[h_idx, f_idx, t_idx] = Tmass[mass_str][h_idx][f_idx][t_idx]
                    end
                end
            end

            count_open_habitat = 0.0
            for i in 1:tmax
                for j in 1:f_bins
                    if T_matrix_for_mass[1, j, i] == 2.0
                        count_open_habitat += 1
                    end
                end
            end

            proportions[mass] = count_open_habitat / (tmax * f_bins)
        else
            println("Aviso: massa $(mass_str) não encontrada.")
        end
    end

    return proportions
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

function fitness_sensibility(Smass, exp_start::Float64, exp_step::Float64, exp_end::Float64)
    fitness_values = mean_fitness(Smass, exp_start, exp_step, exp_end)
    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents
    
    fit_sensibility = Float64[]
    fitness_array = [fitness_values[m] for m in masses]
    
    # Sensibilidade na escala log: ΔF/Δ(log10 M)
    for i in 1:(length(exponents)-1)
        push!(fit_sensibility, (fitness_array[i+1] - fitness_array[i]) / (exponents[i+1] - exponents[i]))
    end
    
    return fit_sensibility
end

function behavior_flexibility(Tmass, exp_start::Float64, exp_step::Float64, exp_end::Float64)
    proportions = prop_open_habitat(Tmass, exp_start, exp_step, exp_end)
    exponents = collect(exp_start:exp_step:exp_end)
    masses = 10 .^ exponents
    
    behavior_flex = Float64[]
    prop_array = [proportions[m] for m in masses]
    
    # Sensibilidade na escala log: ΔZ/Δ(log10 M)
    for i in 1:(length(exponents)-1)
        push!(behavior_flex, (prop_array[i+1] - prop_array[i]) / (exponents[i+1] - exponents[i]))
    end
    
    return behavior_flex
end

function process_fisio_type(fisio_type::String, exp_start::Float64, exp_step::Float64, exp_end::Float64, run_date::String)
    Tmass, Smass = open_files(fisio_type, exp_start, exp_step, exp_end, run_date)
    
    if isempty(Tmass) || isempty(Smass)
        return (Float64[], Float64[], Float64[])
    end
    
    ΔZΔM = behavior_flexibility(Tmass, exp_start, exp_step, exp_end)
    ΔFΔM = fitness_sensibility(Smass, exp_start, exp_step, exp_end)
    exponents = collect(exp_start:exp_step:exp_end)
    
    return (exponents, ΔZΔM, ΔFΔM)
end

# ============================================================
# MAIN
# ============================================================
if length(ARGS) < 7
    println("Uso: julia analyze_SDP.jl <fisio_type_browser> <fisio_type_grazer> <fisio_type_mixed> <exp_start> <exp_end> <exp_step> <run_date>")
    println("Exemplo: julia sensibility_analysis.jl browsers grazers mixed 0 4 0.1 20260908")
    exit(1)
end

fisio_type_browser = ARGS[1]
fisio_type_grazer  = ARGS[2]
fisio_type_mixed   = ARGS[3]
exp_start = parse(Float64, ARGS[4])
exp_end   = parse(Float64, ARGS[5])
exp_step  = parse(Float64, ARGS[6])
run_date  = ARGS[7]

println("\n=== Processando dados ===")

# Processa os três tipos fisiológicos
exponents_b, ΔZΔM_browsers, ΔFΔM_browsers = process_fisio_type(fisio_type_browser, exp_start, exp_step, exp_end, run_date)
exponents_g, ΔZΔM_grazers, ΔFΔM_grazers   = process_fisio_type(fisio_type_grazer, exp_start, exp_step, exp_end, run_date)
exponents_m, ΔZΔM_mixed, ΔFΔM_mixed       = process_fisio_type(fisio_type_mixed, exp_start, exp_step, exp_end, run_date)

# Usa o eixo log10 (remove o primeiro ponto porque diff reduz o tamanho em 1)
masses_log10_b = exponents_b[2:end]
masses_log10_g = exponents_g[2:end]
masses_log10_m = exponents_m[2:end]

println("ΔF/ΔM por massa → Browsers: $(collect(zip(round.(10 .^ masses_log10_b, digits=2), ΔFΔM_browsers))) | Grazers: $(collect(zip(round.(10 .^ masses_log10_g, digits=2), ΔFΔM_grazers))) | Mixed: $(collect(zip(round.(10 .^ masses_log10_m, digits=2), ΔFΔM_mixed)))")

println("\n=== Gerando gráficos combinados ===")
println("Tamanhos: Browsers=$(length(ΔZΔM_browsers)), Grazers=$(length(ΔZΔM_grazers)), Mixed=$(length(ΔZΔM_mixed))")

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

# ============================================================
# Xticks customizados (potências inteiras de 10)
# ============================================================
exp_start_int = Int(round(exp_start))
exp_end_int   = Int(round(exp_end))
xticks_values = collect(exp_start_int:exp_end_int)
xticks_labels = ["10" * join(['\u2070','\u00B9','\u00B2','\u00B3','\u2074','\u2075','\u2076','\u2077','\u2078','\u2079'][d .+ 1] for d in digits(n)) for n in xticks_values]

# ============================================================
# GRÁFICO 1: ΔZ/ΔM (Comportamento) - TODAS AS CURVAS JUNTAS
# ============================================================
# p1 = plot(masses_log10_b, ΔZΔM_browsers,
#           label      = "Browsers",
#           lw         = 5,
#           marker     = :none,
#           color      = browser_color,
#           xlabel     = "Body mass (kg)",
#           ylabel     = "ΔZ / Δ(log₁₀ M)",
#           legend     = :best,
#           size       = (800, 650),
#           grid       = false,
#           xlim       = (exp_start - 0.2, exp_end + 0.2),
#           framestyle = :box,
#           xticks     = (xticks_values, xticks_labels))

# plot!(p1, masses_log10_g, ΔZΔM_grazers,
#       label  = "Grazers",
#       lw     = 5,
#       marker = :none,
#       color  = grazer_color)

# plot!(p1, masses_log10_m, ΔZΔM_mixed,
#       label  = "Mixed Feeders",
#       lw     = 5,
#       marker = :none,
#       color  = mixed_color)

# hline!(p1, [0], color=:black, linestyle=:dash, label="", alpha=0.5, lw=2)

# ============================================================
# GRÁFICO 2: ΔF/ΔM (Fitness) - TODAS AS CURVAS JUNTAS
# ============================================================
p2 = plot(masses_log10_b, ΔFΔM_browsers,
          #label      = "Browsers",
          lw         = 5,
          marker     = :none,
          color      = browser_color,
          xlabel     = "Body mass (kg)",
          ylabel     = "ΔFitness / ΔMass",
          legend     = :none,
          size       = (800, 650),
          grid       = false,
          xlim       = (exp_start - 0.2, exp_end + 0.2),
          framestyle = :box,
          xticks     = (xticks_values, xticks_labels))

plot!(p2, masses_log10_g, ΔFΔM_grazers,
      #label  = "Grazers",
      lw     = 5,
      marker = :none,
      color  = grazer_color)

plot!(p2, masses_log10_m, ΔFΔM_mixed,
      #label  = "Mixed Feeders",
      lw     = 5,
      marker = :none,
      color  = mixed_color)

hline!(p2, [0], color=:black, linestyle=:dash, #label="", 
        alpha=0.5, lw=2)
# vline!(p2, [log10(3.5)],
#     lw =2,
#     linestyle = :dash,
#     color = :red
# )

# vline!(p2, [log10(6.5)],
#     lw =2,
#     linestyle = :dash,
#     color = :red
# )
vspan!(p2, [log10(3.5), log10(10)],
    color = :gray,
    alpha = 0.2,
    label = ""
)
# ============================================================
# Salvando os gráficos
# ============================================================
script_dir = @__DIR__
results_dir = joinpath(script_dir, "results")
isdir(results_dir) || mkpath(results_dir)

#output_behavior = joinpath(results_dir, "ΔZ_combined_$(run_date).svg")
output_fitness  = joinpath(results_dir, "ΔF_combined_$(run_date).svg")


#savefig(p1, output_behavior)
savefig(p2, output_fitness)

#println("✓ Gráfico de comportamento (ΔZ) salvo: $(output_behavior)")
println("✓ Gráfico de fitness (ΔF) salvo: $(output_fitness)")
#println("\n✅ Análises concluídas com sucesso para todos os fisiotipos.")