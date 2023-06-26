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

hub_rename_dict = Dict(
    "HUB1" => "VindØ",
    "HUB2" => "NSWHP",
    "HUB3" => "Bornholm"
)

function select_weeks(df, weeks)
    return df[[w in weeks for w in df[!, "weeks"]], :]
end

function select_mean(df, by, cols)
    return combine(groupby(df, by),
        cols .=> mean .=> cols)
end

function select_mean_weighted(df, by, cols, weights)
    return combine(groupby(df, by),
        cols .=> (c,w) -> mean(c,Weights(w)) .=> cols)
end

function load_electrolyser_results_df(
    model_results,
    hub_rename_dict)
    electrolyser_df = DataFrame(gettable(
        model_results["electrolysers"],
        infer_eltypes=true)...)
    for var in ["DA", "BA_down", "BA_up"]
        electrolyser_df[!, var] = Float64.(
            electrolyser_df[!, var])
    end

    for (key,val) in hub_rename_dict
        electrolyser_df[!,:n] = 
            ifelse.(
                electrolyser_df[!,:n] .== key,
                val,
                electrolyser_df[!,:n])
    end
    
    return electrolyser_df
end

function calculate_electrolyser_period_weights(
    electrolyser_df
)
    electrolyser_df[!, "cap_fac_sum"] = 
        electrolyser_df[!, :DA] .+
        electrolyser_df[!, :BA_down] .+ 
        electrolyser_df[!, :BA_up]
    agg_df = combine(groupby(electrolyser_df, [:scenario, :e]),
        :cap_fac_sum => sum => :cap_fac_sum_agg)
    electrolyser_df = 
        innerjoin(electrolyser_df, agg_df, on = [:scenario, :e])
    electrolyser_df[!, :weights] = 
        electrolyser_df[!, :cap_fac_sum] ./
        electrolyser_df[!, :cap_fac_sum_agg]

    replace!(electrolyser_df.mean_expected_power_price, NaN => 0)
    replace!(electrolyser_df.weights, NaN => 0)
    replace!(electrolyser_df.mean_expected_power_price, "NaN" => 0)
    replace!(electrolyser_df.mean_expected_power_price, 1.0 => 0)

    return electrolyser_df
end


path_results = "results/data/"

timestamp = "2022-06-15-00-00"
filename_aggregated_results = "model_results_"*timestamp*".xlsx"

aggregated_results = openxlsx(path_results*filename_aggregated_results)

model_statistics_df = DataFrame(gettable(
    aggregated_results["model_statistics"],
    infer_eltypes=true))

############### Electrolyser plots ###############
electrolyser_df = 
    load_electrolyser_results_df(
        aggregated_results,
        hub_rename_dict)

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

######## Plotting electrolyser ########
df = copy(electrolyser_df)
weeks = unique(df[!, :weeks])

df[!,:scenario] = df[!,:scenario] .* "_" .* string.(df[!,:no_scenarios])

flt_bz =
    occursin.("TYNDP", df[!, "scenario"]) .&
    occursin.("2030", df[!, "scenario"]) .&
    (df[!, "ntc_scaling_factor"] .== 1.0)

flt_year = 
    occursin.("TYNDP", df[!, "scenario"]) .&
    occursin.("OBZ", df[!, "scenario"]) .&
    (df[!, "ntc_scaling_factor"] .== 1.0)

flt_data = 
    occursin.("OBZ", df[!, "scenario"]) .&
    occursin.("2030", df[!, "scenario"]) .&
    (df[!, "ntc_scaling_factor"] .== 1.0)

flt_ntc = 
    occursin.("OBZ", df[!, "scenario"]) .&
    occursin.("2030", df[!, "scenario"]) .&
    occursin.("TYNDP", df[!, "scenario"])

configurations = [
    Dict(
        "case" => "BZ_comparison",
        "flt" => flt_bz,
        "fname" => "BZ_comparison_2030"),
    Dict(
        "case" => "year_comparison",
        "flt" => flt_year,
        "fname" => "2030_2040_comparison"),
    Dict(
        "case" => "dataset_comparison",
        "flt" => flt_data,
        "fname" => "TYNDP_oE_comparison"),
    Dict(
        "case" => "ntc_comparison",
        "flt" => flt_ntc,
        "fname" => "ntc_comparison")
]

# config = configurations[1]

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

# for week in 1:length(weeks)
#     flt =
#         occursin.("TYNDP", df[!, "scenario"]) .&
#         # occursin.("OBZ", df[!, "scenario"])# .&
#         occursin.("2030", df[!, "scenario"]) .&
#         (df[!,:weeks] .== weeks[week])

#     f = plot_capacity_factors_grouped_from_dataframe(df[flt,:])
#     @show f

#     f |> FileIO.save(
#         "results/plots/cap_fac_grp_$(timestamp)_OBZvsHBZ_$(week).pdf")
# end


########## Calculate electrolyser statistics for selected countries ##########
zones = vcat(["DKW1", "DKE1"], values(hub_rename_dict)...)
flt = [(n in zones) & (y == 2030) & (ntc == 1.0) & (bz .== "OBZ") for (n, y, ntc, bz) in
    zip(electrolyser_df[!, :n],
        electrolyser_df[!, :y],
        electrolyser_df[!, :ntc_scaling_factor],
        electrolyser_df[!, :bz_config])]
df = electrolyser_df[flt,:]
df_mean = select_mean(
    df,
    [:scenario, :n],
    [:DA, :BA_down, :BA_up]
)
df_mean = leftjoin(df_mean, unique(electrolyser_df[!,[:n, :g_max]]), on=:n)
df_mean[!, "cap_fac_sum"] = df_mean[!, :DA] .+ df_mean[!, :BA_down] .+ df_mean[!, :BA_up]
df_mean[!, "prod"] = df_mean[!,"cap_fac_sum"] .* df_mean[!, :g_max] 


##### Select TYNDP
electrolyser_df =
    electrolyser_df[occursin.("TYNDP", electrolyser_df[!, :scenario]), :]

electrolyser_df[electrolyser_df[!,:scenario] .== "OBZ_2030_TYNDP_TYNDP_1.0", :]
# electrolyser_df = electrolyser_df[electrolyser_df[!,:scenario] .== "OBZ_2030_openENTRANCE_Paris_1.0", :]

zones = vcat(["DKW1", "DKE1"], values(hub_rename_dict)...)
flt = [(n in zones) & (y == 2030) & (ntc == 1.0) for (n, y, ntc) in
    zip(electrolyser_df[!, :n], electrolyser_df[!, :y], electrolyser_df[!, :ntc_scaling_factor])]
# flt = [(y == 2030) for y in electrolyser_df[!, :y]]
electrolyser_df = electrolyser_df[flt, :]


############### Manually calculate profit of electrolyser ###############
electrolyser_df[!, "prod"] = (electrolyser_df[!, :DA] .+ electrolyser_df[!, :BA_down] .-
    electrolyser_df[!, :BA_up]) .* electrolyser_df[!, :g_max] .* 8760/6 ./ 1000 # MWh to GWh

electrolyser_df[!, "mc"] = 
    (electrolyser_df[!, :DA] * 10.6875 .+
    electrolyser_df[!, :BA_down] * 10.6875 * 1.2 .-
    electrolyser_df[!, :BA_up] * 10.6875 * 0.2) ./ 
    (electrolyser_df[!, :DA] .+ electrolyser_df[!, :BA_down] .+ electrolyser_df[!, :BA_up]) 

electrolyser_df[!, "revenue"] = electrolyser_df[!, "prod"] * 150

electrolyser_df[!, "profit_total"] =
    electrolyser_df[!, "prod"] .* 
        (150 .- electrolyser_df[!, "mc"] .- electrolyser_df[!, "mean_expected_power_price"])

electrolyser_df[!, "relative_profit"] = 
    electrolyser_df[!, :profit_total] ./ electrolyser_df[!, :g_max]

flt = [(n in values(hub_rename_dict)) & (y == 2030) for (n, y) in
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

df = df[[1,3,2,4], :]

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


# agg_df = combine(groupby(res_df, [:scenario, :e]),
#     :cap_fac_sum => sum => :cap_fac_sum_agg)
# electrolyser_df = 
#     innerjoin(electrolyser_df, agg_df, on = [:scenario, :e])
# electrolyser_df[!, :weights] = 
#     electrolyser_df[!, :cap_fac_sum] ./
#     electrolyser_df[!, :cap_fac_sum_agg]

df = combine(groupby(res_df, [:scenario, :j]),
    :bz_config => first => :bz_config,
    :n => first => :n,
    ["mean_expected_power_price", "weights"] => 
        ((p,w) -> mean(p, Weights(w))) => "mean_expected_power_price",
    :expected_profit => sum => :expected_profit,
    :expected_spill => sum => :expected_spill,
)
# select!(df, Not([:scenario, :j]))

# electrolyser_df =
#     electrolyser_df[occursin.("TYNDP", electrolyser_df[!, :scenario]), :]

# df = select_mean(
#     electrolyser_df,
#     [:scenario, :e],
#     [:DA, :BA_down, :BA_up]
# )

unique(df[!, :weeks])

df = copy(electrolyser_df)
# copy(select_weeks(electrolyser_df, ["2:2"]))

weeks = unique(df[!, :weeks])
for week in 1:length(weeks)
    flt = occursin.("openENTRANCE", df[!, "scenario"]) .&
        # occursin.("OBZ", df[!, "scenario"])# .&
        occursin.("2030", df[!, "scenario"]) .&
        (df[!,:weeks] .== weeks[week])

    f = plot_capacity_factors_grouped_from_dataframe(df[flt,:])
    @show f

    f |> FileIO.save(
        "results/plots/capacity_factors_grouped_$(timestamp)_OBZvsHBZ_oE_$(week).pdf")
end

############### Timeseries plots ###############
scenarios = [
    "OBZ_2030_TYNDP_TYNDP_1.0",
    "HBZ_2030_TYNDP_TYNDP_1.0",
    ]
results_ts_dict = Dict(
    sc => CSV.read(path_results*"output_"*sc*"_"*timestamp*".csv",
    DataFrame, delim=",") for sc in scenarios
)


sum(ts[!, "price_DA_HUB1"] .* ts[!, "L_HUB1_ptg_offshore"] 
./ sum(ts[!, "L_HUB1_ptg_offshore"]))


########## Price duration curve ##########
s_list = [
    "OBZ_2030_TYNDP_TYNDP_1.0",
    "HBZ_2030_TYNDP_TYNDP_1.0",
    ]
f = plot_price_duration_curve_from_dataframe(
    "HUB1", Dict(s => results_ts_dict[s] for s in s_list))
f |> FileIO.save(
    "results/plots/price_duration_curve_$timestamp.pdf")

f = plot_electrolyser_duration_curve_from_dataframe(
    "HUB1_ptg_offshore", Dict(s => results_ts_dict[s] for s in s_list))
f |> FileIO.save(
    "results/plots/electrolyser_duration_curve_$timestamp.pdf")

Interconnectors_df = DataFrame(gettable(
    openxlsx(path*"ntc-v18-yw-2022_05_16.xlsx")["ntc"])...)

I = [
    # "DKW1_HUB1",
    "HUB1_DKW1",
    # "HUB3_DKE1"
    # "DKE1_HUB3"
    ]
s_list = [
    "OBZ_2030_TYNDP_TYNDP_1.0",
    "HBZ_2030_TYNDP_TYNDP_1.0"
    ]
f = plot_interconnector_flow_from_dataframe(
    I,
    get_ntcs(I, Interconnectors_df),
    Dict(s => results_ts_dict[s] for s in s_list))
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
#2755 + 72
k = 2335
f = plot_system_balance_from_dataframe(
    scenario_res,
    k:k+72,#150:150+96,
    results_ts_dict["OBZ_2030_TYNDP_TYNDP_1.0"])
@show f
f |> FileIO.save(
    "results/plots/system_balance_$timestamp.pdf")

node_df = DataFrame(gettable(
    aggregated_results["countries"],
    infer_eltypes=true)...)
flt = [(y == 2030) & (occursin("TYNDP", sc) & (bz .== "OBZ")) 
    for (y,sc,bz) in zip(node_df[!, :y], node_df[!, :scenario], node_df[!, :bz_config])]
node_df = node_df[flt, :]

node_df[(node_df[!, :scenario] .== "OBZ_2030_TYNDP_TYNDP_1.0") .&
    (node_df[!,:n] .=="DKE1") , "mean_power_price_DA"]

#300:400
f = plot_RES_generation("DELU_onshore_wind", 310:382, ES.Ω, ES)
f |> FileIO.save(
    "results/plots/wind_gen.pdf")

#4058:4102
f = plot_RES_generation("DELU_solarpv", 4008:4080, ES.Ω, ES)
f |> FileIO.save(
    "results/plots/pv_gen.pdf")


############ Number scenarios evaluation ############
sc_obj_val_df = DataFrame(gettable(
    aggregated_results["scenario_obj_vals"],
    infer_eltypes=true)...)

sc_obj_val_df[!, :scenario] = 
    sc_obj_val_df[!, :scenario] .* "_" .* 
    string.(sc_obj_val_df[!, :no_scnearios])

cols = vcat("scenario", "weeks", "s" .* string.(collect(1:10)))
sc_obj_val_df[!,cols]

using KernelDensity, Random, Statistics

f = plot()
week = "3:6"
for sc in unique(sc_obj_val_df[!, :scenario])
    
    x = convert(Array{Float64,1}, 
    collect(sc_obj_val_df[
        (sc_obj_val_df[!, :scenario] .== sc) .&
        (sc_obj_val_df[!, :weeks] .== week), cols[3:end]][1,:]))
    println("sc: $sc", " mean: $(mean(x))", " var: $(std(x))")
    y = kde(x)
    plot!(f, 
        # range(minimum(x),step=10000,stop=maximum(x)),
        range(4.90e9, step=100000, stop=5.35e9),
        z->pdf(y, z),
        label=sc)
end
@show f