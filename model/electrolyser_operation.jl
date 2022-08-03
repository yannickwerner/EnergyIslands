# Set file directory as working directory
cd(dirname(@__FILE__))
using Pkg
Pkg.activate("")

using JuMP, Gurobi
using CSV, DataFrames
using XLSX: openxlsx, gettable, writetable
using Random: seed!
using StatsBase: sample, mean, Weights
using Statistics: quantile
using Dates
using DataStructures: counter
using Impute: interp
using Plots, VegaLite, FileIO

include("components.jl")
include("data.jl")
include("energy_system.jl")
include("opt_model.jl")
include("post_processing.jl")
include("plotting.jl")

# Implement week/hours information/time windows properly

# week_list = [1:8760]

week_list = [1:10]

# h = 1460
# week_list = [
#     1:h,
#     h+1:h*2,
#     h*2+1:h*3,
#     h*3+1:h*4,
#     h*4+1:h*5,
#     h*5+1:h*6
# ]

config_dict = Dict(
    :case_study => "29032022",
    :year_pointforecast => 2018, # one in 2015:2019
    :years_scenarios => 2015:2019, #2016:2019
    :countries =>
        ["IE00", "DKE1", "DKW1", "DELU", "DKKF", "DEKF",
         "BE00", "FR00", "NL00", "HUB1", "HUB2", "HUB3",
         "NOM1", "NON1", "NOS0", "UK00", "SE01", "SE02",
         "SE03", "SE04", "PL00", "EE00", "LT00", "LV00"],
    :μ => 0.2,
    :seed => seed!(03042022),
    :hydrogen_autarky_rate => 0.6,
    :quantile_ramping => 0.999,
    :wind_uncertainty_scaling_factor => 1.5,
    :hydrogen_production_driver => "price", # One of "price" or "demand
    # :co2_price_pathway => "Ambitious", # One of "TYNDP", "BAU", "Ambitious", "Paris"
    # :data_set_capacities => "TYNDP", # One of ["TYNDP", "openENTRANCE"]
    # :ntc_scaling_factor => 1.0, 
    # :uncertainty_horizon => "long-termfuture" # "near-termfuture", "long-termfuture"
    # :bidding_zone_config => "OBZ" # One of "OBZ", "HBZ"
)

scenarios = [
    ("OBZ", 2030, "TYNDP", "TYNDP", 1.0),
    # ("OBZ", 2030, "openENTRANCE", "Paris", 1.0),
    # ("OBZ", 2040, "TYNDP", "TYNDP", 1.0),
    # ("OBZ", 2040, "openENTRANCE", "Paris", 1.0),
    # ("HBZ", 2030, "TYNDP", "TYNDP", 1.0),
    # ("HBZ", 2030, "openENTRANCE", "Paris", 1.0),
    # ("HBZ", 2040, "TYNDP", "TYNDP", 1.0),
    # ("HBZ", 2040, "openENTRANCE", "Paris", 1.0),
    # ("OBZ", 2030, "PAUL", "BAU", 1.0),
    # ("OBZ", 2040, "PAUL", "BAU", 1.0),
    # ("HBZ", 2030, "PAUL", "BAU", 1.0),
    # ("HBZ", 2040, "PAUL", "BAU", 1.0),
    # ("OBZ", 2030, "TYNDP", "TYNDP", 0.6),
    # ("OBZ", 2040, "TYNDP", "TYNDP", 0.6),
    # ("OBZ", 2030, "TYNDP", "TYNDP", 0.8),
    # ("OBZ", 2030, "TYNDP", "TYNDP", 1.2),
    # ("OBZ", 2040, "TYNDP", "TYNDP", 0.8),
    # ("OBZ", 2040, "TYNDP", "TYNDP", 1.2),
    ]

results = Dict()

# Preparations for postprocessing
timestamp = Dates.format(Dates.now(), "yyyy-mm-dd-HH-MM")
countries_data_extracted =
    ["HUB1", "HUB2", "HUB3", "DKW1", "DKE1", "DELU"]
scenarios_data_extraction = sample(
    config_dict[:seed],
    ["s$i" for i in 1:length(config_dict[:years_scenarios])],
    1,
    replace=false)
global electrolyser_df = DataFrame()
global renewables_df = DataFrame()
global node_df = DataFrame()
global model_statistics_df = DataFrame()
global scenario_obj_val_df = DataFrame()

for sc in scenarios
    config_dict[:scenario] = join(sc, "_")
    config_dict[:bidding_zone_config] = sc[1]
    config_dict[:year_model] = sc[2]
    config_dict[:data_set_capacities] = sc[3]
    config_dict[:co2_price_pathway] = sc[4]
    config_dict[:ntc_scaling_factor] = sc[5]

    global optimal_values_ts_df = DataFrame()
    for weeks in week_list    
        # timesteps = Int.(round.(reduce(vcat, [collect((k-1)*168+1:k*168) 
        #     for k in weeks])))
        timesteps = collect(weeks)
        config_dict[:T] = timesteps
        global ES = EnergySystem(config_dict...)
        load_input_data!(ES)
        ES.T = 1:length(ES.T)
        start_time = Dates.now()
        run_opt_model!(ES; print_model=false);
        end_time = Dates.now()

        ##### Postprocessing #####
        calculate_descriptive_statistics!(ES);
        ES.total_run_time = string(Dates.canonicalize(
            Dates.CompoundPeriod(end_time-start_time)))
        results[sc[1]*"_"*string(sc[2])*"_"*string(sc[3])] = ES
        ##### Extract optimal solutions
        e_df, r_df, n_df, m_s_df, obj_val_df = 
            write_descriptive_statistics(ES, weeks)
        global electrolyser_df = vcat(electrolyser_df, e_df)
        global renewables_df = vcat(renewables_df, r_df)
        global node_df = vcat(node_df, n_df)
        global model_statistics_df = vcat(model_statistics_df, m_s_df)
        global scenario_obj_val_df = vcat(scenario_obj_val_df, obj_val_df)
        df = extract_optimal_values(
            scenarios_data_extraction,
            countries_data_extracted,
            ES)
        df[!, :t] = timesteps
        global optimal_values_ts_df = vcat(
            optimal_values_ts_df, df
        )
    end
    CSV.write(
        "results/data/output_$(
            join(sc[1:5], "_"))_$timestamp.csv",
        optimal_values_ts_df)
end

writetable(
    "results/data/model_results_$timestamp.xlsx",
    "model_statistics" => model_statistics_df,
    "electrolysers" => electrolyser_df,
    "renewables" => renewables_df,
    "countries" => node_df,
    "scenario_obj_vals" => scenario_obj_val_df
    )


