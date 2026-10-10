#!/usr/bin/env julia
using JSON3
using Plots
using StatsPlots
using StatsBase

# ==============================================================================
# ARGUMENTOS DA LINHA DE COMANDO
# ==============================================================================

if length(ARGS) < 4
    println("Uso: julia plot_SDP.jl <fisio_type> <data> <mass> <patch>")
    println("Exemplo: julia plot_SDP.jl grazers 20250114 1500 1")
    println("  fisio_type: tipo fisiológico (grazers, browsers, etc)")
    println("  data: data no formato YYYYMMDD")
    println("  mass: massa corporal (g)")
    println("  patch: número do patch (1, 2, etc)")
    exit(1)
end

fisio_type = ARGS[1]
date_str = ARGS[2]
mass = parse(Float64, ARGS[3])
h = parse(Int, ARGS[4])

# ==============================================================================
# CARREGAR DADOS
# ==============================================================================

script_dir = @__DIR__
data_dir = joinpath(script_dir, "data")

# Constrói o nome do arquivo baseado nos argumentos
input_file = joinpath(data_dir, "SDP_$(fisio_type)_$(date_str).json")

if !isfile(input_file)
    println("Erro: Arquivo '$(input_file)' não encontrado.")
    println("Procurando em: $(data_dir)")
    exit(1)
end

println("Carregando dados de: $(input_file)")
data = JSON3.read(input_file, Dict)

# ==============================================================================
# CONVERTER DADOS PARA ARRAYS
# ==============================================================================

# Converte as estruturas aninhadas de volta para arrays 3D
function nested_to_array3d(nested_vec)
    h_size = length(nested_vec)
    y_size = length(nested_vec[1])
    t_size = length(nested_vec[1][1])
    
    arr = zeros(Float64, h_size, y_size, t_size)
    for i in 1:h_size
        for j in 1:y_size
            for k in 1:t_size
                arr[i, j, k] = nested_vec[i][j][k]
            end
        end
    end
    return arr
end

# Carrega S e T
S_dict = Dict()
T_dict = Dict()

for (mass_str, nested_data) in data["S"]
    S_dict[parse(Float64, mass_str)] = nested_to_array3d(nested_data)
end

for (mass_str, nested_data) in data["T"]
    T_dict[parse(Float64, mass_str)] = nested_to_array3d(nested_data)
end

println("Massas disponíveis: ", sort(collect(keys(S_dict))))

# ==============================================================================
# VERIFICAR SE A MASSA EXISTE
# ==============================================================================

mass_rounded = round(mass; digits=2)
if !haskey(S_dict, mass_rounded)
    println("Erro: Massa $(mass_rounded) não encontrada nos dados.")
    println("Massas disponíveis: ", sort(collect(keys(S_dict))))
    exit(1)
end

# ==============================================================================
# EXTRAIR DADOS PARA PLOTAGEM
# ==============================================================================

S_array = S_dict[mass_rounded]
T_array = T_dict[mass_rounded]

# Dimensões
tmax = size(S_array, 3)
f_bins = size(S_array, 2)

println("\nDimensões dos dados:")
println("  Patches: ", size(S_array, 1))
println("  Fat bins: ", f_bins)
println("  Tempo: ", tmax)

# Seleciona o slice para o patch específico
S_slice = S_array[h, :, :]
T_slice = T_array[h, :, :]

# ==============================================================================
# PLOTAGEM
# ==============================================================================

# Plot 1: Survival (S)
plot_S = StatsPlots.heatmap(
    1:tmax, 1:f_bins, S_slice,
    xlabel="Time",
    ylabel="Fat State (bin)",
    title="Survival Probability - Mass=$(mass_rounded)g, Patch=$h",
    colorbar_title="S(h,y,t)",
    c=:inferno,
    size=(700, 500)
)

# Plot 2: Strategy (T)
# Identifica valores únicos (patches) na matriz
unique_patches = sort(unique(T_slice))
n_patches = length(unique_patches)

# Cria colormap discreto
colors = cgrad(:viridis, n_patches, categorical=true)

plot_T = StatsPlots.heatmap(
    1:tmax, 1:f_bins, T_slice,
    xlabel="Time",
    ylabel="Fat State (bin)",
    title="Optimal Strategy - Mass=$(mass_rounded)g, Patch=$h",
    colorbar_title="Patch Choice",
    c=colors,
    clims=(minimum(unique_patches) - 0.5, maximum(unique_patches) + 0.5),
    size=(700, 500)
)

# Plot combinado
plot_combined = plot(
    plot_S, plot_T,
    layout=(1, 2),
    size=(1400, 500),
    plot_title="SDP Results - $(fisio_type)"
)

display(plot_combined)

# ==============================================================================
# SALVAR PLOTS (OPCIONAL)
# ==============================================================================

output_dir = joinpath(data_dir, "plots")
if !isdir(output_dir)
    mkdir(output_dir)
end

savefig(plot_S, joinpath(output_dir, "S_$(fisio_type)_mass$(mass_rounded)_patch$(h).png"))
savefig(plot_T, joinpath(output_dir, "T_$(fisio_type)_mass$(mass_rounded)_patch$(h).png"))
savefig(plot_combined, joinpath(output_dir, "combined_$(fisio_type)_mass$(mass_rounded)_patch$(h).png"))

println("\n✓ Plots salvos em: $(output_dir)")
println("  - Survival (S)")
println("  - Strategy (T)")
println("  - Combined")

# ==============================================================================
# ANÁLISE ADICIONAL (OPCIONAL)
# ==============================================================================

# Estatísticas básicas
println("\n=== ESTATÍSTICAS ===")
println("Survival (S):")
println("  Min: ", minimum(S_slice))
println("  Max: ", maximum(S_slice))
println("  Mean: ", round(sum(S_slice) / length(S_slice); digits=4))

println("\nStrategy (T):")
println("  Patches únicos escolhidos: ", sort(unique(T_slice)))
println("  Patch mais comum: ", mode(T_slice[T_slice .> 0]))  # Ignora zeros