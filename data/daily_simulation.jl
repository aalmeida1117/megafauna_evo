using Distributions
using JSON3
using Dates
using LinearAlgebra
using Base.Threads 

include("allometric_functions.jl")
using .AllometricFunctions
##com reaction plane

#hello 
function dailysim(mass, patch, fisio_type, mu)
    ### 1 = closed habitat
    ### 2 = open habitat
    velocity = find_velocity(mass) # Velocidade do forrageador em m/s
    beta = bite_size_allo(mass)

    if patch == 1 
        velocity = 0.5*velocity # Ajuste para habitat fechado
    end

    if fisio_type == "grazers" && patch == 1
        beta = 0.5*beta # Ajuste para habitat fechado
    end
    
    kcal_final = 0.0
    #active_hours = foraging_hours(mass)
    tmax_bout, _ = foragingtime(mass) .* (60 * 60)
    _, n = indperarea(mass)
    chewrate = chew_allo(mass) # g/s
    tchewgram = 1 / chewrate # s/g
    tchew = tchewgram * beta # s/g * g/bite = s/bite
    width = reactionwidth(mass)
    height = reactionheight(mass)
    b0_bmr = 3.4 # Taxa metabólica basal constante (watts kg^-0.75)
    b0_fmr = 6.8 # Taxa metabólica de campo constante (watts kg^-0.75)
    basal_mr = b0_bmr * (mass ^ 0.75) # Taxa metabólica basal para kg
    field_mr = b0_fmr * (mass ^ 0.75) # Taxa metabólica ativa para kg
    upper = 18.0 #maior e menores valores de x na normal
    lower = 0.0 #maior e menores valores de x na normal

    t_travel = 0.0
    t_chew = 0.0
    bites = 0.0
    GUT = 0.0

    total_kcal = 0.0
    cost_ws = 0.0
    number_of_successes = 0.0
    t = 0.0
    gut = 0.0

    max_gut = gut_volume_g(mass)

    if patch == 2 && height>2
            height = 2# Ajuste para habitat aberto
    end

    if fisio_type == "mixed"
        #res_traits p/ grazers
        res_traits_g=  [mu=mu alpha= 4.5 mean_edensity = mean_energy_closed(mass) sd_edensity = 2.8 zeta= 2.2 d_pred = 0.0002;
                        mu= mu alpha = 9 mean_edensity = mean_energy_open(mass) sd_edensity = 1 zeta = 1 d_pred = 0.0002;];

        mu_g = res_traits_g[patch,1]
        alpha_g= res_traits_g[patch, 2]
        mean_edensity_g = res_traits_g[patch, 3]
        sd_edensity_g = res_traits_g[patch, 4]
        zeta_g = res_traits_g[patch, 5]

        #res_traits p/ browser
        res_traits_b=  [mu=mu alpha= 9 mean_edensity = mean_energy_closed(mass) sd_edensity = 2.8 zeta= 1 d_pred = 0.0002;
                        mu=mu alpha = 4.5 mean_edensity = mean_energy_open(mass) sd_edensity = 1 zeta = 2.2 d_pred = 0.0002;];               
        
        mu_b = res_traits_b[patch, 1]
        alpha_b = res_traits_b[patch, 2]
        mean_edensity_b = res_traits_b[patch, 3]
        sd_edensity_b = res_traits_b[patch, 4]
        zeta_b = res_traits_b[patch, 5]
        d_pred = res_traits_b[patch,6]
        #edensity = Normal(mean_edensity, sd_edensity)
        #Reaction plane area
        
        if patch == 2 && height>2
            height = 2# Ajuste para habitat aberto
        end
                        #Competitor density (per reaction plane area)
        na = n/(width*height);
                        
                        # Bite density field
        m_res_b = mu_b * (1 / beta) #* height #* width
        m_res_g = mu_g *(1/beta)

                        #Rescaling bite density by Competitor Load
        mprime_b = m_res_b / na #m para browser
        mprime_g = m_res_g / na #m para grazers

        alphaprime_b = alpha_b * na^(zeta_b - 2) #alpha para browsers
        alphaprime_g = alpha_g * na^(zeta_g - 2) #alpha para browsers

        
        edensity_b = Truncated(Normal(mean_edensity_b, sd_edensity_b), lower, upper)
        edensity_g = Truncated(Normal(mean_edensity_g, sd_edensity_g), lower, upper)

        # Inicialização de variáveis escalares


        gammadist_b = Gamma(alphaprime_b, mprime_b / alphaprime_b)
        gammadist_g = Gamma(alphaprime_g, mprime_g / alphaprime_g)

        while t < tmax_bout
            distance_to_resource_b = rand(Exponential(1.0 / rand(gammadist_b)))
            distance_to_resource_g = rand(Exponential(1.0 / rand(gammadist_g)))

            deltat_b = distance_to_resource_b / velocity
            cost_ws_b = field_mr * deltat_b
            cost_kcal_b = cost_ws_b*0.0002

            deltat_g = distance_to_resource_g / velocity
            cost_ws_g = field_mr * deltat_g
            cost_kcal_g = cost_ws_g*0.0002
            # Atualiza tempo e custos
            kcal_food_b = rand(edensity_b)
            kcal_food_g = rand(edensity_g)

            if (kcal_food_b-cost_kcal_b)>(kcal_food_g-cost_kcal_g)

                deltat = deltat_b
                t += deltat
                t_travel += deltat
                cost_ws += field_mr * deltat

                if (tmax_bout > (t + deltat)) && (max_gut > gut)
                    number_of_successes += 1
                    
                    # tb = total_biomass(mass, res_traits, patch)
                    # if Int64(round(tb/ bite_size_allo(mass))) < 0
                    #     println("negative total biomass")
                    # end
                    #try
                    bites += 1
                    #catch e 
                    #    display("deu erro", e)
                    #end
                    
                    t += tchew
                    t_chew += tchew 
                    gut += beta
                    cost_ws += field_mr * tchew + d_pred
                    GUT = gut
                    kcal_food = kcal_food_b
                    kcal_final += beta * kcal_food
                else
                    break
                end
            else
                deltat = deltat_g
                t += deltat
                t_travel += deltat
                cost_ws += field_mr * deltat

                if (tmax_bout > (t + deltat)) && (max_gut > gut)
                    number_of_successes += 1
                    
                    # tb = total_biomass(mass, res_traits, patch)
                    # if Int64(round(tb/ bite_size_allo(mass))) < 0
                    #     println("negative total biomass")
                    # end
                    #try
                    bites += 1
                    #catch e 
                    #    display("deu erro", e)
                    #end
                    
                    t += tchew
                    t_chew += tchew 
                    gut += beta
                    cost_ws += field_mr * tchew + d_pred
                    GUT = gut
                    kcal_food = kcal_food_g
                    kcal_final += beta * kcal_food
                else
                    break
                end
            end
        end

    else
        if fisio_type == "grazers" 
            res_traits= [mu=mu alpha= 4.5 mean_edensity = mean_energy_closed(mass) sd_edensity = 2.8 zeta= 2.2 d_pred = 0.0002;
                        mu= mu alpha = 9 mean_edensity = mean_energy_open(mass) sd_edensity = 1 zeta = 1 d_pred = 0.0002;];
        
        
        elseif fisio_type == "browsers" 
            res_traits= [mu=mu alpha= 9 mean_edensity = mean_energy_closed(mass) sd_edensity = 2.8 zeta= 1 d_pred = 0.0002;
                        mu=mu alpha = 4.5 mean_edensity = mean_energy_open(mass) sd_edensity = 1 zeta = 2.2 d_pred = 0.0002;];
        
        end

        # Mover as definições para antes do uso
        mu = res_traits[patch, 1]
        alpha = res_traits[patch, 2]
        mean_edensity = res_traits[patch, 3]
        sd_edensity = res_traits[patch, 4]
        zeta = res_traits[patch, 5]
        d_pred = res_traits[patch,6]
        #edensity = Normal(mean_edensity, sd_edensity)
        #Reaction plane area
        width = reactionwidth(mass)
        height = reactionheight(mass)

        if patch == 2 && height>2
            height = 2# Ajuste para habitat aberto
        end
                        #Competitor density (per reaction plane area)
        na = n/(width*height);
                        
                        # Bite density field
        m_res = mu * (1 / beta) #* height #* width

                        #Rescaling bite density by Competitor Load
        mprime = m_res / na

        alphaprime = alpha * na^(zeta - 2)
        # WATTS/KG
        edensity = Truncated(Normal(mean_edensity, sd_edensity), lower, upper)
        gammadist = Gamma(alphaprime, mprime / alphaprime)
        

        while t < tmax_bout
            distance_to_resource = rand(Exponential(1.0 / rand(gammadist)))

            # Atualiza tempo e custos
            deltat = distance_to_resource / velocity
            t += deltat
            t_travel += deltat
            cost_ws += field_mr * deltat

            if (tmax_bout > (t + deltat)) && (max_gut > gut)
                number_of_successes += 1
                
                # tb = total_biomass(mass, res_traits, patch)
                # if Int64(round(tb/ bite_size_allo(mass))) < 0
                #     println("negative total biomass")
                # end
                #try
                bites += 1
                #catch e 
                #    display("deu erro", e)
                #end
                
                t += tchew
                t_chew += tchew 
                gut += beta
                cost_ws += field_mr * tchew + d_pred
                GUT = gut
                kcal_food = rand(edensity)
                kcal_final += beta * kcal_food
            else
                break
            end
        end
    end
    total_t = 60 * 60 * 24
    rest_t = total_t - t
    cost_ws += basal_mr * rest_t
    cost_whr = cost_ws / 3600 # Convertendo watt/s para watt/hr
    cost_kcal = cost_whr * 0.86 # Convertendo watt para kcal
    total_kcal = kcal_final

    return total_kcal, number_of_successes, cost_kcal
end

# ========== PARÂMETROS VIA TERMINAL ==========
if length(ARGS) < 5
    println("Uso: julia daily_simulation.jl exp_start exp_end exp_step n_configs fisio_type")
    println("Ex: julia daily_simulation.jl 1 3 0.25 20000 grazers")
    exit(1)
end

exp_start = parse(Float64, ARGS[1])   # expoente inicial (ex: 1 → 10^1)
exp_end   = parse(Float64, ARGS[2])   # expoente final   (ex: 3 → 10^3)
exp_step  = parse(Float64, ARGS[3])   # passo entre expoentes
n_configs = parse(Int, ARGS[4])       # número de simulações
fisio_type = ARGS[5]                  # "grazers", "browsers" ou "mixed"

nbins = 50
mass_dict = Dict()
habitats = [1, 2]  # 1 = fechado, 2 = aberto

# Gera as massas como potências de 10
exponents = collect(exp_start:exp_step:exp_end)
masses = 10 .^ collect(exp_start:exp_step:exp_end)
mu_values = range(0.0001, 0.001, length = 30)

# Garante que a pasta existe antes das threads escreverem
isdir("data") || mkpath("data")


# println("\n=== PARÂMETROS DA SIMULAÇÃO ===")
# println("Expoente inicial: ", exp_start)
# println("Expoente final: ", exp_end)
# println("Passo: ", exp_step)
# println("Expoentes: ", exponents)
# println("Número de expoentes: ", length(exponents))
# println("\nMassas geradas (kg):")
# for (i, m) in enumerate(masses)
#     println("  $i. 10^$(exponents[i]) = $(round(m, digits=2)) kg")
# end
println("\nNúmero total de massas: ", length(masses))
println("Tipo fisio: ", fisio_type)
# println("Configurações por massa/habitat: ", n_configs)
# println("="^40)

# ========== FUNÇÕES AUXILIARES ==========

function calculate_histogram(data, nbins)
    min_val = minimum(data)
    max_val = maximum(data)
    bin_edges = range(min_val, max_val, length=nbins+1)
    bin_counts = zeros(Int, nbins)

    for val in data
        for i = 1:nbins
            if val >= bin_edges[i] && val < bin_edges[i+1]
                bin_counts[i] += 1
                break
            end
        end
    end

    if data[end] == max_val
        bin_counts[end] += 1
    end

    return bin_edges, bin_counts
end

function find_bin_index(value, bins)
    nbins = length(bins) - 1
    for i = 1:nbins
        if value >= bins[i] && value < bins[i+1]
            return i
        end
    end
    return nbins
end

# ========== LOOP DE SIMULAÇÃO (escala log10) ==========

println("\n=== INICIANDO SIMULAÇÕES ===")

@threads for mu in mu_values
    println("mu=$(round(mu, sigdigits=4)) rodando na thread $(Threads.threadid())")
    mu_str = string(round(mu; sigdigits = 4))
    mass_dict = Dict()
    for (idx, mass) in enumerate(masses)
        println("\n" * "="^50)
        println("MASSA $(idx)/$(length(masses)): $(round(mass; digits=4)) kg")
        println("Expoente: 10^$(exponents[idx])")
        println("="^50)
        
        # Converte massa para string para usar como chave (evita problemas de precisão float)
        mass_key = string(round(mass; digits=2))
        mass_dict[mass_key] = Dict()

        for patch in habitats
            # patch_name = patch == 1 ? "fechado" : "aberto"
            # println("\n  → Habitat $(patch_name) (patch=$patch)")
            
            mass_dict[mass_key][patch] = Dict()
            gains = Float64[]
            costs = Float64[]

            for config in 1:n_configs
                if config % 5000 == 0
                    # print("    Progresso: $(config)/$(n_configs)\r")
                    flush(stdout)
                end
                
                try
                    total_mcal, _, cost_mcal = dailysim(mass, patch, fisio_type,mu)
                    push!(gains, total_mcal)
                    push!(costs, cost_mcal)
                catch e
                    println("\n    ❌ ERRO na config $config:")
                    println("       ", e)
                    rethrow(e)
                end
            end
            # println("    Progresso: $(n_configs)/$(n_configs) ✓")
            
            #println("    Calculando histogramas...")
            gains_bin_edges, gains_bin_counts = calculate_histogram(gains, nbins)
            costs_bin_edges, costs_bin_counts = calculate_histogram(costs, nbins)

            gains_probabilities = gains_bin_counts / sum(gains_bin_counts)
            costs_probabilities = costs_bin_counts / sum(costs_bin_counts)

            gains_bin_midpoints = (gains_bin_edges[1:end-1] + gains_bin_edges[2:end]) / 2
            costs_bin_midpoints = (costs_bin_edges[1:end-1] + costs_bin_edges[2:end]) / 2

            gain_cost_pairs = [(gains[i], costs[i]) for i in 1:length(gains)]
            gain_bins = range(minimum(gains), maximum(gains), length=nbins)
            cost_bins = range(minimum(costs), maximum(costs), length=nbins)

            joint_histogram = zeros(Int, nbins, nbins)

            for pair in gain_cost_pairs
                gain_bin_index = find_bin_index(pair[1], gain_bins)
                cost_bin_index = find_bin_index(pair[2], cost_bins)
                joint_histogram[gain_bin_index, cost_bin_index] += 1
            end

            joint_probabilities = joint_histogram / sum(joint_histogram)

            mass_dict[mass_key][patch]["gains"] = collect(gain_bins)
            mass_dict[mass_key][patch]["costs"] = collect(cost_bins)
            mass_dict[mass_key][patch]["prob"]  = joint_probabilities
            
            #println("    ✓ Habitat $(patch_name) concluído")
        end
        
        println("\n✓ Massa $(mass_key) kg processada completamente")
    end
    output_file = joinpath("data", "$(fisio_type)_mu$(mu_str)_dict.json")
    JSON3.write(output_file, mass_dict)
    println("✓ mu=$mu_str salvo em $output_file")
end

# ========== SALVAR SAÍDA ==========
# isdir("data") || mkpath("data")
# output_file = joinpath("data", "$(fisio_type)_mu$(mu_str)_dict.json")
# JSON3.write(output_file, mass_dict)

#println("\nSimulação finalizada! Resultados salvos em $(output_file).")


## como rodar:
#julia --threads 80 daily_simulation_mu.jl 1 3 0.25 20000 grazers
# ou, se tiver menos núcleos:
#julia --threads auto daily_simulation_mu.jl 1 3 0.25 20000 grazers
