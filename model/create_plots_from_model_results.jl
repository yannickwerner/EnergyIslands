# Set file directory as working directory
cd(dirname(@__FILE__))
using Pkg
Pkg.activate("")

using CSV, DataFrames
using XLSX: openxlsx, gettable, writetable
using Plots#; pyplot()
using StatsPlots
using Plots.PlotMeasures 
using VegaLite, FileIO
using StatsBase: mean, Weights

include("plotting.jl")
include("data.jl")

# bz_rename_dict = Dict(
#     "HUB1" => "DEI",
#     "HUB2" => "NSWHP",
#     "HUB3" => "Bornholm"
# )

bz_rename_dict = Dict(
    "HUB1" => "DEI",
    "HUB2" => "NSWHP",
    "HUB3" => "Bornholm",
    "DKW1" => "DK1",
    "DKE1" => "DK2",
    "NL00" => "NL",
    "DELU" => "DE",
    "BE00" => "BE",
    "UK00" => "UK",
    "PL00" => "PL",
    "SE01" => "SE1",
    "SE02" => "SE2",
    "SE03" => "SE3",
    "SE04" => "SE4",
    "NON1" => "NON",
    "NOS0" => "NOS",
    "NOM1" => "NOM",
    "FR00" => "FR"
)

path_results = "results/data/"

# timestamp = "2022-06-15-00-00"
timestamp = "2023-05-10-00-00"
filename_aggregated_results = "model_results_"*timestamp*".xlsx"

aggregated_results = openxlsx(path_results*filename_aggregated_results)

model_statistics_df = DataFrame(gettable(
    aggregated_results["model_statistics"],
    infer_eltypes=true))

############### Electrolyser plots ###############
electrolyser_df = 
    load_electrolyser_results_df(
        aggregated_results,
        bz_rename_dict)

##### Add electrolyser capacities
path = "data/29032022/"
e_df = DataFrame(gettable(
    openxlsx(path*"inst_cap-v15-yw-2022_06_08.xlsx")["g_max"]))
e_df[!, "name"] = [n*"_"*g for (n,g) in 
    zip(e_df[!, "n"], 
    e_df[!, "g"])]
e_df = e_df[e_df[!, :y] .== 2030, :]
electrolyser_df[!, "g_max"] .= 0.0
for row in 1:nrow(electrolyser_df)
    try
        electrolyser_df[row, :g_max] = 
            e_df[e_df[!, :name] .== electrolyser_df[row, :e], :g_max][1]
    catch Error
        continue
    end
end

electrolyser_df = calculate_electrolyser_period_weights(electrolyser_df);

############## Scenario number comparison ##############
filename_aggregated_results = "model_results_2023-05-12-00-00_no_scenario_comp.xlsx"
model_statistics_df = DataFrame(gettable(
    openxlsx(path_results*filename_aggregated_results)["electrolysers"],
    infer_eltypes=true))
df_mean = select_mean(
    model_statistics_df,
    [:no_scenarios, :n],
    [:DA, :BA_down, :BA_up]
)

flt = [(n in ["HUB1", "HUB3"]) for  n in df_mean[!, :n]]
df_mean[flt, :]


######## Plotting electrolyser ########
df = copy(electrolyser_df)
weeks = unique(df[!, :weeks])

# Add number of uncertainty scenarios to case string 
df[!,:scenario] = df[!,:scenario] .* "_" .* string.(df[!,:no_scenarios])


configurations = define_filters(df)
# config = configurations[4]

for config in configurations

    df_mean = select_mean(
        df[config["flt"],:],
        [:scenario, :n],
        [:DA, :BA_down, :BA_up]
    )

    f = plot_capacity_factors_grouped_from_dataframe(
        df_mean, config["case"])
    @show f

    f |> FileIO.save(
        "results/plots/cap_fac_grp_$(timestamp)_$(config["fname"]).pdf")
end


# Comparison of different weeks
# f = compare_weeks(df, weeks, timestamp)
# @show f

# ########## Calculate electrolyser statistics for selected countries ##########
# zones = vcat(["DKW1", "DKE1"], values(bz_rename_dict)...)

zones = values(bz_rename_dict)
flt = [(occursin("TYNDP_base", sc)) & (n in zones) & (y == 2030) & (bz .== "OBZ") for (sc, n, y, bz) in
    zip(electrolyser_df[!, :scenario],
        electrolyser_df[!, :n],
        electrolyser_df[!, :y],
        # electrolyser_df[!, :ntc_scaling_factor],
        electrolyser_df[!, :bz_config])]
        #flt = [(n in zones) & (y == 2030) & (ntc == 1.0) & (bz .== "OBZ") for (n, y, ntc, bz)
df = electrolyser_df[flt,:]
df_mean = select_mean(
    df,
    [:scenario, :n],
    [:DA, :BA_down, :BA_up]
)

df_mean = leftjoin(df_mean, unique(electrolyser_df[!,[:n, :g_max]]), on=:n)
df_mean[!, "cap_fac_sum"] = df_mean[!, :DA] .+ df_mean[!, :BA_down] .+ df_mean[!, :BA_up]
df_mean[!, "prod"] = df_mean[!,"cap_fac_sum"] .* df_mean[!, :g_max] 

df_mean[!, :share_BA] = 
    (df_mean[!, :BA_down] .- df_mean[!, :BA_up]) ./ df_mean[!, :DA]

##### Select TYNDP
electrolyser_df =
    electrolyser_df[
        occursin.("openENTRANCE", electrolyser_df[!, :scenario]) .& #TYNDP
        occursin.("base", electrolyser_df[!, :scenario]), :]

# electrolyser_df = electrolyser_df[electrolyser_df[!,:scenario] .== "OBZ_2030_TYNDP_TYNDP_base_1.0", :]
# electrolyser_df = electrolyser_df[electrolyser_df[!,:scenario] .== "OBZ_2030_openENTRANCE_Paris_1.0", :]

zones = vcat(["DKW1", "DKE1"], values(bz_rename_dict)...)
# zones = unique(electrolyser_df[!,:n])
flt = [(n in zones) & (y == 2030) & (ntc == 1.0) for (n, y, ntc) in
    zip(electrolyser_df[!, :n], electrolyser_df[!, :y], electrolyser_df[!, :ntc_scaling_factor])]
# flt = [(y == 2030) for y in electrolyser_df[!, :y]]
electrolyser_df = electrolyser_df[flt, :]


############### Manually calculate profit of electrolyser ###############
electrolyser_df[!, "prod"] = (electrolyser_df[!, :DA] .+ electrolyser_df[!, :BA_down] .-
    electrolyser_df[!, :BA_up]) .* electrolyser_df[!, :g_max] .* 8760/6 ./ 1000 # MWh to GWh

electrolyser_df[!, "mc"] = 
    (electrolyser_df[!, :DA] * 10.6875 .+
    electrolyser_df[!, :BA_down] * 10.6875 * 1.2 .+
    electrolyser_df[!, :BA_up] * 10.6875 * 0.8) ./ 
    (electrolyser_df[!, :DA] .+ electrolyser_df[!, :BA_down] .+ electrolyser_df[!, :BA_up]) 

electrolyser_df[!, "revenue"] = electrolyser_df[!, "prod"] * 150

electrolyser_df[!, "profit_total"] =
    electrolyser_df[!, "prod"] .* 
        (150 .- electrolyser_df[!, "mc"] .- electrolyser_df[!, "mean_expected_power_price"])

electrolyser_df[!, "relative_profit"] = 
    electrolyser_df[!, :profit_total] ./ electrolyser_df[!, :g_max]

flt = [(n in values(bz_rename_dict)) & (y == 2030) for (n, y) in
    zip(electrolyser_df[!, :n], electrolyser_df[!, :y])]

df = electrolyser_df[flt, 
    ["scenario", "e", "weeks", "bz_config", "n", "mean_expected_power_price", 
     "profit_total", "prod", "weights", "mc"]]
df[!, "profit_total"] = df[!, "profit_total"]./1e6

# use weighted averae for power price
df = combine(groupby(df, [:scenario, :e]),
    :bz_config => first => :bz_config,
    :n => first => :n,
    ["mean_expected_power_price", "weights"] => 
        ((p,w) -> mean(p, Weights(w))) => "mean_expected_power_price",
    ["mc", "weights"] => 
        ((mc,w) -> mean(mc, Weights(w))) => "mc_weighted",
    :profit_total => sum => :profit_total,
    :prod => sum => :prod,
)


df[!, ["mean_expected_power_price", "profit_total", "prod", "mc_weighted"]] = round.(
    df[!, ["mean_expected_power_price", "profit_total", "prod", "mc_weighted"]],
        digits=2)

df = df[!, 
    ["bz_config", "n", "prod", "mc_weighted", "mean_expected_power_price", "profit_total"]]

# df = df[[1,3,2,4], :]

df[!, "profit_total"] = 
    df[!, "prod"] .* (
        150 .- df[!, "mc_weighted"]/0.665 .- df[!,"mean_expected_power_price"]/0.665)

########## Business case manual calculation ##########
using Latexify
latexify(df; env=:table, latex=false)

println("OBZ:", 
sum(electrolyser_df[electrolyser_df[!,:bz_config] .== "OBZ", "prod"]))

println("HBZ:", 
sum(electrolyser_df[electrolyser_df[!,:bz_config] .== "HBZ", "prod"]))

# electrolyser_df[electrolyser_df[!, "scenario"] .== "OBZ_2030_TYNDP", :]

##### WindFarms
zones = ["HUB1", "HUB3"]#, "DKW1", "DKE1"]
res_df = DataFrame(gettable(
    aggregated_results["renewables"],
    infer_eltypes=true)...)
flt = [(n in zones) & (y == 2030) & occursin("TYNDP", sc) & (ntc == 1.0)
    for (n, y, sc, ntc) in
    zip(res_df[!, :n], res_df[!, :y], res_df[!, :scenario], res_df[!, :ntc_scaling_factor])]
res_df = res_df[flt,:]
res_df[!, "expected_profit"] = res_df[!, "expected_profit"] ./ 1e6
res_df[!, "expected_spill"] = res_df[!, "expected_spill"] ./ 1e3

names(res_df)


df = combine(groupby(res_df, [:scenario, :j]),
    :bz_config => first => :bz_config,
    :n => first => :n,
    ["mean_expected_power_price", "weights"] => 
        ((p,w) -> mean(p, Weights(w))) => "mean_expected_power_price",
    :expected_profit => sum => :expected_profit,
    :expected_spill => sum => :expected_spill,
)

############### Timeseries plots ###############
s_list = [
    "OBZ_2030_TYNDP_TYNDP_base_1.0",
    "HBZ_2030_TYNDP_TYNDP_base_1.0",
    ]
results_ts_dict = Dict(
    sc => CSV.read(path_results*"output_"*sc*"_"*timestamp*".csv",
    DataFrame, delim=",") for sc in s_list
)

########## Price duration curve ##########
f = plot_price_duration_curve_from_dataframe(
    "HUB1", Dict(s => results_ts_dict[s] for s in s_list))
f |> FileIO.save(
    "results/plots/price_duration_curve_$timestamp.pdf")

f = plot_electrolyser_duration_curve_from_dataframe(
    "HUB1_ptg_offshore", Dict(s => results_ts_dict[s] for s in s_list))
f |> FileIO.save(
    "results/plots/electrolyser_duration_curve_$timestamp.pdf")


########## Interconnector flows ##########
Interconnectors_df = DataFrame(gettable(
    openxlsx(path*"ntc-v18-yw-2022_05_16.xlsx")["ntc"]))

I = [# "DKW1_HUB1",# "HUB3_DKE1"# "DKE1_HUB3"
    "HUB1_DKW1"]
f = plot_interconnector_flow_from_dataframe(
    I, get_ntcs(I, Interconnectors_df), Dict(s => results_ts_dict[s] for s in s_list))
@show f
f |> FileIO.save(
    "results/plots/interconnector_flow_shore_to_EI_$timestamp.pdf")

i = "HUB1_DKW1"
f = plot_price_during_congestion_from_dataframe(
    i,
    2030,
    Dict(s => results_ts_dict[s] for s in s_list),
    get_ntcs([i], Interconnectors_df);
    T=100:172)
@show f

f |> FileIO.save(
    "results/plots/price_during_congestion_$timestamp.pdf")

########## ##########
i = "HUB1_DKW1"
result_dict = Dict(s => results_ts_dict[s] for s in s_list)
ntc_dict = get_ntcs([i], Interconnectors_df)
year = 2030
T = 1:72
for sc in keys(result_dict)
    expected_trade = result_dict[sc][!, "F_$(i)"] .+
        result_dict[sc][!, "F_adj_$(i)_avg"]

    hours_congestion = expected_trade .>= ntc_dict[string(year)][i]
    if isnothing(T)
        T = 1:length(hours_congestion[hours_congestion .== 1])
    end

    for c in rsplit(i, "_")

        if ((c == "HUB1") & (rsplit(sc,"_")[1] == "HBZ"))
            continue
        end

        println("$sc, $c, $(mean(
            result_dict[sc][!, "price_DA_$c"][hours_congestion][T]))")
    end
end


names_df = names(results_ts_dict[collect(keys(results_ts_dict))[1]])
scenario_res = rsplit(names_df[occursin.(
    "system_balancing_conventionals", names_df)][1], "_")[end]
# k = 2335
j = 241
f = plot_system_balance_from_dataframe(
    scenario_res,
    j:j+72,
    # k:k+72,
    results_ts_dict["OBZ_2030_TYNDP_TYNDP_base_1.0"])
@show f
f |> FileIO.save(
    "results/plots/system_balance_$timestamp.pdf")


node_df = DataFrame(gettable(
    aggregated_results["countries"],
    infer_eltypes=true)...)
flt = [(y == 2030) & (occursin("TYNDP", sc) & (bz .== "OBZ")) 
    for (y,sc,bz) in zip(node_df[!, :y], node_df[!, :scenario], node_df[!, :bz_config])]
node_df = node_df[flt, :]

node_df[(node_df[!, :scenario] .== "OBZ_2030_TYNDP_TYNDP_base_1.0") .&
    (node_df[!,:n] .=="DKE1") , "mean_power_price_DA"]

#300:400
f = plot_RES_generation("DELU_onshore_wind", 310:382, ES.Ω, ES)
f |> FileIO.save(
    "results/plots/wind_gen.pdf")

#4058:4102
f = plot_RES_generation("DELU_solarpv", 4008:4080, ES.Ω, ES)
f |> FileIO.save(
    "results/plots/pv_gen.pdf")


# ############ Number scenarios evaluation ############
# sc_obj_val_df = DataFrame(gettable(
#     aggregated_results["scenario_obj_vals"],
#     infer_eltypes=true)...)

# sc_obj_val_df[!, :scenario] = 
#     sc_obj_val_df[!, :scenario] .* "_" .* 
#     string.(sc_obj_val_df[!, :no_scnearios])

# cols = vcat("scenario", "weeks", "s" .* string.(collect(1:10)))
# sc_obj_val_df[!,cols]

# using KernelDensity, Random, Statistics

# f = plot()
# week = "3:6"
# for sc in unique(sc_obj_val_df[!, :scenario])
    
#     x = convert(Array{Float64,1}, 
#     collect(sc_obj_val_df[
#         (sc_obj_val_df[!, :scenario] .== sc) .&
#         (sc_obj_val_df[!, :weeks] .== week), cols[3:end]][1,:]))
#     println("sc: $sc", " mean: $(mean(x))", " var: $(std(x))")
#     y = kde(x)
#     plot!(f, 
#         # range(minimum(x),step=10000,stop=maximum(x)),
#         range(4.90e9, step=100000, stop=5.35e9),
#         z->pdf(y, z),
#         label=sc)
# end
# @show f