using XLSX
using StatsPlots
using DataFrames

#ler o doc em xlsx
df = DataFrame(XLSX.readtable("pnas201907998_s5_movfsn.xlsx", 1))
#verificar nome das colunas
names(df)

