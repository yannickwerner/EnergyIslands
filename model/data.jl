function bound_vector(vec, ub, lb)
    """
    Returns a vector that is lower and upper bounded by
    given bounds by replacing values outside of those bounds.

    Parameters
    ----------
    vec : Array
        Vector of floats

    ub : float
        Upper bound
    
    lb : float
        Lower bound

    Returns
    -------
    vec : Array
        Vector with replaced bounds
    """
    vec = ifelse.(vec .>= ub, ub, vec)
    vec = ifelse.(vec .<= lb, lb, vec)
    return vec
end

function get_ntcs(I, I_df)
    """
    Retrieves net-transfer capacities (NTCs) for given 
    interconnectors from the input data for 2030 and 2040.

    Parameters
    ----------
    I : Array
        List of interconnectors

    I_df : DataFrame
        Raw input data containing all NTCs

    Returns
    -------
    ntc_dict : Dict
        Contains NTCs for given interconnectors
    """
    ntc_dict = Dict()
    for year in [2030, 2040]
        ntc_dict_year = Dict()
        for i in I
            flt = [(
                (rsplit(i, "_")[1] == n) &&
                (rsplit(i, "_")[2] == m) &&
                (y == year)) ?
                true : false for (n,m,y) in eachrow(I_df[!, [:n, :m, :y]])]
            ntc_dict_year[i] = I_df[flt, :ntc][1]
        end
        ntc_dict[string(year)] = ntc_dict_year
    end
    return ntc_dict
end

function extract_ts_data_pv(ES, df, country, tech, g_max)
    """
    Retrieves and reformats pv generation time series data
    for a given bidding zone. Uses data from different
    years as scenarios to create a probabilistic forecast 
    for a single year. 

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    df : DataFrame
        Raw input data for pv time series for all countries
        and time steps

    country : str
        Bidding zone name

    tech : str
        Technology for which time series data is retrieved
        (pv only)

    g_max : float
        Installed capacity of pv in bidding zone used to Normalize
        time series data

    Returns
    -------
    point_forecast : Array
        Contains bidding zone specific point forecast data for each time steps
    
    scenario_output : Dict
        Contains probabilistic pv forecast
    """
    if isnothing(ES.years_scenarios)
        scenario_list = collect(1980:2019)
        scenario_output = Dict("s" * string(y-1979) =>
            g_max .* df[
                (df[!, :y] .== y) .&
                (df[!, :technology] .== tech), country][ES.T]
                    for y in scenario_list)
    else
        scenario_list = ES.years_scenarios
        scenario_dict_omega = Dict(scenario_list[y] => 
            "s" * string(y) for y in 1:length(scenario_list))
        scenario_output = Dict(scenario_dict_omega[y] =>
            g_max .* df[
                (df[!, :y] .== y) .&
                (df[!, :technology] .== tech), country][ES.T]
                    for y in scenario_list)
    end

    flt = [((t in ES.T) && (y in scenario_list)) ?
        true : false for (t,y) in eachrow(df[!, [:t, :y]])]
    point_forecast = g_max .* combine(groupby(df[flt, :], :t), 
        country => mean => country)[!,country]

    return point_forecast, scenario_output
end


function extract_ts_data_wind(ES, df, country, tech, g_max)
    """
    Retrieves and reformats wind generation time series data
    for a given bidding zone. Uses forecast errors based on
    forecast and real-time realizations from different years
    and uses those as scenarios for a single year.

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    df : DataFrame
        Raw input data for wind time series for all countries
        and time steps

    country : str
        Bidding zone name

    tech : str
        Technology for which time series data is retrieved
        and reformatted (wind only)

    g_max : float
        Installed capacity of wind (onshore/offshore) in
        bidding zone used to normalize time series data

    Returns
    -------
    point_forecast : Array
        Contains bidding zone specific point forecast data for each time steps
    
    scenario_output : Dict
        Contains probabilistic wind forecast
    """
    years = 1980:2020
    sc_dict = Dict(years[y] => 
        "s" * string(y) for y in 1:(2020-1980))
    scenario_list = ES.years_scenarios

    col = country * "_" * tech

    pf_normalized = df[
        (df[!, :y] .== ES.year_pointforecast) .&
        (df[!, :technology] .== tech), country][ES.T]

    scenario_output = Dict("s"*string(y) =>
        bound_vector(
            g_max .* (pf_normalized .+
            ES.forecast_error_scenarios[col][
                sc_dict[scenario_list[y]]]), g_max, 0)
        for y in 1:length(scenario_list))

    point_forecast = vec(mean(Array(DataFrame(scenario_output)), dims=2))

    return point_forecast, scenario_output
end


function extract_ts_data_ror(ES, df, country, tech, g_max)
    """
    Retrieves and reformats run-of-river generation time
    series data for a given bidding zone. Uses data from 
    different years as scenarios to create a probabilistic
    forecast for a single year. Since data is only available
    from 2014--2019, other years are randomly sampled from
    those years given.

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    df : DataFrame
        Raw input data for pv time series for all countries
        and time steps

    country : str
        Bidding zone name

    tech : str
        Technology for which time series data is retrieved
        (pv only)

    g_max : float
        Installed capacity of pv in bidding zone used to Normalize
        time series data

    Returns
    -------
    point_forecast : Array
        Contains bidding zone specific point forecast data for each time steps
    
    scenario_output : Dict
        Contains probabilistic pv forecast
    """

    years_scenario = ES.hydro_sample_years["years_scenarios"]

    scenario_output = Dict("s" * string(i) =>
        g_max .* df[
            (df[!, :y] .== years_scenario[i]) .&
            (df[!, :technology] .== tech), country][ES.T]
                for i in 1:length(years_scenario))

    flt = [((t in ES.T) && (y in ES.hydro_sample_years["years_unique"])) ?
        true : false for (t,y) in eachrow(df[!, [:t, :y]])]
    point_forecast = g_max .* combine(groupby(df[flt, :], :t), 
        [country, "weights"] => (
            (c,w) -> mean(c, Weights(w))) => country)[!, country]

    return point_forecast, scenario_output
end


function extract_tech_subset(df, techs, ES)
    flt = [((g in techs) && (y == ES.year_model) && (n in ES.countries)) ?
        true : false for (g,y,n) in eachrow(df[!, [:g, :y, :n]])]
    return df[flt, [:n, :g, :g_max]]
end


function calculate_generation_limit(ES, df, country)
    """
    Calculates total generation limit for hydro reservoir production
    based on historical time series data.

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    df : DataFrame
        Contains historical production data for hydro reservoirs 
        for all countries and time steps
    
    country : str
        Bidding zone name

    Returns
    -------
    scenario_output : Dict
        Contains probabilistic forecast of total available
        power production of hydro reservoir power
    """
    years_scenario = ES.hydro_sample_years["years_scenarios"]

    scenario_output = Dict()
    for i in 1:length(years_scenario)
        flt = [((t in ES.T) && (y == years_scenario[i])) ?
            true : false for (t,y) in eachrow(df[!, [:t, :y]])]
        scenario_output["s" * string(i)] = sum(df[flt, country])
    end
            
    return scenario_output
end


########## Caculate sample years for hydro ##########
# If the year selected for
# the point forecast is not available, i.e. in 2015--2019, then
# the year is randomly drawn from the available ones.
function select_hydro_sample_years(RunOfRiver_ts, ES)
    """
    Calculate probabilistic forecast for run-of-river technology
    and samples years for the probabilistic modeling of available
    hydro reservoir power production based on historical data.
    Ensures that hydro reservoir and run-of-river timeseries are consistent.

    Parameters
    ----------
    RunOfRiver_ts : DataFrame
        Contains historical production data for run-of-river technology
        for all countries and time steps

    ES : EnergySystem
        Contains information about energy system modeled

    Returns
    -------
    RunOfRiver_ts : DataFrame
        Contains weights for each year to calculate expected
        power production from run-of-river
    
    hydro_sample_years : Dict
        Contains the randomly sampled years for the probabilistic
        forecast of hydro power production
    """
    ### Available years in data ###
    available_years = collect(2015:2019)
    ### Scenario years ###
    # Determine sample years for hydro based on available data
    if isnothing(ES.years_scenarios)
        years_scenario_hydro = sample(ES.seed, available_years, 2014-1979)
        # Assure that the last 5 years coincide with the real ones
        append!(years_scenario_hydro, available_years)
    else
        # Replace years that are not available in the data
        # with randomly sampled years.
        scenario_list_years = ES.years_scenarios
        unavailable_years = [~(y in available_years) ? true : false
            for y in scenario_list_years]
        sample_years = sample(ES.seed, available_years, sum(unavailable_years))
        years_scenario_hydro = scenario_list_years .* (.~unavailable_years)

        # overwrite given years without timeseries data
        j = 1
        for i in 1:length(years_scenario_hydro)
            if years_scenario_hydro[i] == 0
                years_scenario_hydro[i] = sample_years[j]
                j += 1
            end
        end
    end

    # Years to exchange scenario data
    year_counter = counter(years_scenario_hydro)
    unique_sample_years_sorted = sort(collect(keys(year_counter)))
    hydro_sample_years_distribution = 1/sum(values(year_counter)) .* [
        year_counter[y] for y in unique_sample_years_sorted
    ]

    # Construct a Dataframe with year weights to determine
    # weighted mean production
    weight_df = DataFrame(
        y = unique_sample_years_sorted,
        weights = hydro_sample_years_distribution
    )
    RunOfRiver_ts = innerjoin(RunOfRiver_ts, weight_df, on=:y)

    hydro_sample_years = Dict{Any,Any}(
        "years_scenarios" => years_scenario_hydro,
        "years_unique" => unique_sample_years_sorted
    )

    return RunOfRiver_ts, hydro_sample_years
end
    

function calculate_forecast_error_scenarios(
    df,
    bidding_zone,
    T,
    scaling_factor)
    """
    Calculate probabilistic forecast errors based on historical
    day-ahead forecasts and real-time realizations of wind power.

    Parameters
    ----------
    df : DataFrame
        Contains forecasts and real-time realizations for all
        countries, technologies, and timesteps

    bidding_zone : str
        Bidding zone name

    T : Array
        List containing relevant time steps

    scaling_factor : float
        Scales wind forecast by 1-scaling_factor percent

    Returns
    -------
    scenario_dict : Dict
        Scaled, probablistic forecast error for given
        technology and country
    """

    cols = names(df)

    # Assign bidding zones to countries if a country has multiple bidding zones
    if bidding_zone in ["DKE1", "DKW1", "HUB1", "DKKF", "DEKF"]
        country = "DK"
    elseif bidding_zone == "HUB2"
        country = "UK00"
    elseif occursin("NO", bidding_zone)
        country = "NO"
    elseif ((occursin("SE", bidding_zone)) | (bidding_zone == "HUB3"))
        country = "SE"
    elseif bidding_zone == "LV00"
        country = "LT00"
    else
        country = bidding_zone
    end

    scenario_dict = Dict()
    cols_country = cols[ifelse.(occursin.(country, cols), true, false)]
    years = 1980:2020

    norm_vec = scaling_factor .* ones(size(df, 1))

    if (any(occursin.("onshore", cols_country)) & 
           any(occursin.("offshore", cols_country)))

        horizon_onshore = ifelse(
            country*"_wind_onshore_near-termfuture" in names(df),
            "_wind_onshore_near-termfuture",
            "_wind_onshore_long-termfuture")
        horizon_offshore = ifelse(
            country*"_wind_offshore_near-termfuture" in names(df),
            "_wind_offshore_near-termfuture",
            "_wind_offshore_long-termfuture")
        
        forecast_error =
            (df[!,country*horizon_onshore] .-
                df[!,country*"_wind_onshore_current"]) .* norm_vec
        scenario_dict[bidding_zone*"_onshore_wind"] = 
            Dict("s"*string(ω) => forecast_error[df[!,:y] .== years[ω]][T]
                for ω in 1:(2020-1980))
        forecast_error =
            (df[!,country*horizon_offshore] .-
                df[!,country*"_wind_offshore_current"]) .* norm_vec
        scenario_dict[bidding_zone*"_offshore_wind"] = 
            Dict("s"*string(ω) => forecast_error[df[!,:y] .== years[ω]][T]
                for ω in 1:(2020-1980))
    else
        horizon = ifelse(
            country*"_wind_national_near-termfuture" in names(df),
            "_wind_national_near-termfuture",
            "_wind_national_long-termfuture")

        forecast_error =
            (df[!,country*horizon] .-
                df[!,country*"_wind_national_current"]) .* norm_vec 

            scenario_dict[bidding_zone*"_onshore_wind"] = 
                Dict("s"*string(ω) => forecast_error[df[!,:y] .== years[ω]][T]
                    for ω in 1:(2020-1980))
            scenario_dict[bidding_zone*"_offshore_wind"] = 
                scenario_dict[bidding_zone*"_onshore_wind"]
        
    end

    return scenario_dict
end


function load_entsoe_generation_data(country, path, fuels)
    """
    Loads and preprocesses generation data from ENTSO-E.

    Removes hour with NAs for annual time shift
    (summer/winter) and interpolates where data contains NAs.

    Parameters
    ----------
    country : str
        Bidding zone name

    path : str
        Relative path to the input file

    fuels : Dict
        Mapping of fuel to column names

    Returns
    -------
    df : DataFrame
        The manipulated input DataFrame
    """

    df = CSV.read(joinpath(
        path, "entsoe_generation_$(country)_26012020.csv"), DataFrame, delim=",")
    select!(df, cat(Symbol.(values(fuels)), [:MTU], dims=1))

    if size(df, 1) == 8760 * 4 + 1 * 4
        times_remove = [
            "26.03.2017 02:00 - 26.03.2017 02:15 (CET)",
            "26.03.2017 02:15 - 26.03.2017 02:30 (CET)",
            "26.03.2017 02:30 - 26.03.2017 02:45 (CET)",
            "26.03.2017 02:45 - 26.03.2017 03:00 (CET)",
            "26.03.2017 02:00 - 26.03.2017 02:15 (CET/CEST)",
            "26.03.2017 02:15 - 26.03.2017 02:30 (CET/CEST)",
            "26.03.2017 02:30 - 26.03.2017 02:45 (CET/CEST)",
            "26.03.2017 02:45 - 26.03.2017 03:00 (CET/CEST)"
            ]

        filter!(row -> !(row.MTU in times_remove),  df)
        df[!, :group] = repeat(1:8760, inner=4)

    elseif size(df, 1) == 8760 * 2 + 1 * 2
        times_remove = [
            "26.03.2017 02:00 - 26.03.2017 02:30 (CET)",
            "26.03.2017 02:30 - 26.03.2017 03:00 (CET)",
            "26.03.2017 02:00 - 26.03.2017 02:30 (CET/CEST)",
            "26.03.2017 02:30 - 26.03.2017 03:00 (CET/CEST)",
            ]

        filter!(row -> !(row.MTU in times_remove),  df)
        df[!, :group] = repeat(1:8760, inner=2)

    elseif size(df, 1) == 8760 + 1
        times_remove = [
            "26.03.2017 02:00 - 26.03.2017 03:00 (CET)",
            "26.03.2017 02:00 - 26.03.2017 03:00 (CET/CEST)"
            ]
        filter!(row -> !(row.MTU in times_remove),  df)
        df[!, :group] = 1:8760

    else
        ArgumentError("Length of DataFrame is inappropriate.")
    end

    select!(df, Not(:MTU))

    df = ifelse.(ismissing.(df), "N/A", df)
    df = ifelse.(df .== "n/e", "N/A", df)
    df = ifelse.(df .== "N/A", missing, df)

    for col in names(df)
        if typeof(df[!, col]) == Vector{Int64}
            nothing
        elseif typeof(df[!, col]) == Vector{Union{Missing, Int64}}
            nothing
        else
            df[!, col] = passmissing(parse).(Float64, df[!, col])
        end
    end

    df = combine(groupby(df, :group),
        Symbol.(values(fuels)) .=> mean .=> Symbol.(values(fuels)));

    df = ifelse.(ismissing.(df), 0, df)
    select!(df, Not(:group));
    
    return df
end


function calculate_historical_technical_constraints(countries, q, path)
    """
    Calculates fuel-specific ramping rates based on historical aggregated
    country data based on  ENTSO-E for the year 2017.

    Parameters
    ----------
    countries : Array
        List of countries (bidding zones)

    q : float
        Percentile for CDF of ramping rates  

    path : str
        Relative path to the input file

    Returns
    -------
    gradient_data : Dict
        Fuel and country-specific ramping rates
    """

    fuels = Dict(
        "biomass" => "Biomass  - Actual Aggregated [MW]",
        "lignite" => "Fossil Brown coal/Lignite  - Actual Aggregated [MW]",
        "gas" => "Fossil Gas  - Actual Aggregated [MW]",
        "hardcoal" => "Fossil Hard coal  - Actual Aggregated [MW]",
        "oil" => "Fossil Oil  - Actual Aggregated [MW]",
        "nuclear" => "Nuclear  - Actual Aggregated [MW]",
        "reservoir" => "Hydro Water Reservoir  - Actual Aggregated [MW]",
        "peat" => "Fossil Peat  - Actual Aggregated [MW]",
        "phs_charge" => "Hydro Pumped Storage  - Actual Consumption [MW]",
        "phs_discharge" => "Hydro Pumped Storage  - Actual Aggregated [MW]",
    )

    bidding_zones_to_drop = ["HUB1", "HUB2", "HUB3", "DKKF", "DEKF"]
    filter!(c -> ~(c in bidding_zones_to_drop), countries)

    g = Dict(country => load_entsoe_generation_data(
        country, path, fuels) for country in countries)
    gradient_data = Dict()
    for c in countries
        country_data_gradient = Dict()
        for k in keys(fuels)
            if any(g[c][!,fuels[k]] .!= 0)
                country_data_gradient[k] = 
                    quantile(abs.(diff(g[c][!,fuels[k]][g[c][!,fuels[k]] .!= 0])), q)/
                        maximum(g[c][!,fuels[k]][g[c][!,fuels[k]] .!= 0])
                if k == "biomass"
                    country_data_gradient[k] = 0.8
                end
            else
                nothing
            end
        end
        if (("phs_discharge" in keys(country_data_gradient)) &
            ~("phs_charge" in keys(country_data_gradient)))
            country_data_gradient["phs_charge"] =
                country_data_gradient["phs_discharge"]
        end
        gradient_data[c] = country_data_gradient
    end
    
    return gradient_data
end


function calculate_chp_min_load(countries, year, path, T; offset=0.3)
    """
    Calculates minimum load for chp plants based on historical
    residential heat demand data from 2013.
    Applies and offset to the time series to account for households
    that are not supplied by chps.

    Parameters
    ----------
    countries : Array
        List of countries (bidding zones)
    year : int
        Basis year for historical data
    path : str
        Relative path to the input file
    T : Array
        TImesteps to calculate load for
    offset : float
        Offset used to shift time series upward.
        Default value 0.3

    Returns
    -------
    min_load_ts_dict : Dict
        CHP min load for given countries
    """

    heat_load_ts = CSV.read(
        joinpath(path, "when2heat_filtered_flt.csv"),
        DataFrame,
        delim=";")

    dateformat =Dates. DateFormat("yyyy-mm-ddTHH:MM:SSZ")
    heat_load_ts[!,:time] = Dates.DateTime.(
        heat_load_ts[!,:utc_timestamp], dateformat)
    heat_load_ts[!,:y] = Dates.year.(heat_load_ts[!,:time])

    year = 2013

    bidding_zones_to_drop = ["HUB1", "HUB2", "HUB3", "DKKF", "DEKF"]
    filter!(c -> ~(c in bidding_zones_to_drop), countries)

    min_load_ts_dict = Dict()
    for country in countries

        if country == "DELU"
            bz = "DE"
        elseif country == "CH00"
            bz = "AT"
        elseif country == "UK00"
            bz = "IE"
        elseif country in ["DKW1", "DKE1"]
            bz = "DK"
        elseif country in ["SE01", "SE02", "SE03", "SE04"]
            bz = "SE"
        elseif country in ["NOS0", "NOM1", "NON1"]  
            bz = "SE"
        else
            bz = country[1:2]
        end

        min_load_ts =
            heat_load_ts[heat_load_ts[!, :y] .== year,
                bz*"_heat_demand_total"]./
            maximum(heat_load_ts[heat_load_ts[!, :y] .== year,
                bz*"_heat_demand_total"]) .+
            offset
    
        min_load_ts_dict[country] =
            ifelse.(min_load_ts .>= 1, 1, min_load_ts)[T]
    end

    return min_load_ts_dict
end


function load_input_data!(ES)
    """
    Loads all necessary input data and creates all objects
    needed for running the optimization model.

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    Returns
    -------
    ES : EnergySystem
        Contains all objects needed for the opimization model
    """

    ############ File path ############
    path = "./data/" * ES.case_study * "/"

    ########## Define Sets ##########
    Sets_df = DataFrame(gettable(
        openxlsx(path*"sets-v07-al-2022_05_10.xlsx")["sets"])...)
    set_countries = filter(!ismissing, Sets_df[!, :n])
    if isnothing(ES.countries)
        ES.countries = set_countries
    else
        for c in ES.countries
            if ~(c in set_countries)
                throw(ArgumentError("$(c) is not in the predefined country set."))
            end
        end
    end
    ES.N = ES.countries
    techs_dispatchable = filter(!ismissing, Sets_df[!, :i])
    techs_res = filter(!ismissing, Sets_df[!, :j])
    techs_ptg = filter(!ismissing, Sets_df[!, :e])
    techs_storage = filter(!ismissing, Sets_df[!, :s])

    ######################### Load input data #########################
    ########## Installed Capacities ##########
    ########## Storage data ##########
    if ES.data_set_capacities == "TYNDP"
        Capacity_df = DataFrame(gettable(
            openxlsx(path*"inst_cap-v15-yw-2022_06_08.xlsx")["g_max"])...)
        Storages_df = DataFrame(gettable(
            openxlsx(path*"storage-v05-yw-2022_04_13.xlsx")["storage"])...)
        Interconnectors_df = DataFrame(gettable(
            openxlsx(path*"ntc-v18-yw-2022_05_16.xlsx")["ntc"])...)
    elseif ES.data_set_capacities == "openENTRANCE"
        filenmame_openentrance = "openentrance-v12-yw-2022_06_08.xlsx"
        Capacity_df = DataFrame(gettable(
            openxlsx(path*filenmame_openentrance)[
                "dir_trans_inst_cap"])...)
        Storages_df = DataFrame(gettable(
            openxlsx(path*filenmame_openentrance)[
                "dir_trans_stor"])...)
        Interconnectors_df = DataFrame(gettable(
            openxlsx(path*"ntc-v18-yw-2022_05_16.xlsx")["ntc"])...)
    elseif ES.data_set_capacities == "PAUL"
        filenmame_paul = "paul_case-v01-2022_05_20.xlsx"
        Capacity_df = DataFrame(gettable(
            openxlsx(path*filenmame_paul)[
                "inst_cap_paul"])...)
        Interconnectors_df = DataFrame(gettable(
            openxlsx(path*filenmame_paul)[
                "ntc_paul"])...)
        Storages_df = DataFrame(gettable(
            openxlsx(path*"storage-v05-yw-2022_04_13.xlsx")["storage"])...)
    else
        throw(ArgumentError("Data set for installed capacities is
            $(ES.data_set_capacities) but must be one of
            'openENTRANCE', 'TYNDP', or 'PAUL'."))
    end

    ########## Technology parameters ##########
    xls_tech = openxlsx(path*"technologies-v20-yw-2022_05_12.xlsx")
    Parameters_i = DataFrame(gettable(xls_tech["i"])...)
    Parameters_j = DataFrame(gettable(xls_tech["j"])...)
    Parameters_e = DataFrame(gettable(xls_tech["e"])...)
    Parameters_s = DataFrame(gettable(xls_tech["s"])...)    
    ########## Commodity prices ##########
    FuelPrice_df = DataFrame(gettable(xls_tech["p_fuel"])...)
    CO2Price_df = DataFrame(gettable(xls_tech["p_co2"])...)
    ########## Timeseries data ####################
    Demand_EL_ts = CSV.read(joinpath(
        path, "load-v12-yw-2022_04_02_flt.csv"), DataFrame, delim=";")
    PV_ts = CSV.read(joinpath(
        path, "solarpv-v04-yw-2022_05_11_flt.csv"), DataFrame, delim=";")
    WindOnshore_ts = CSV.read(joinpath(
        path, "wind_onshore-v06-yw-2022_05_11_flt.csv"), DataFrame, delim=";")
    WindOffshore_ts = CSV.read(joinpath(
        path, "wind_offshore-v07-yw-2022_05_11_flt.csv"), DataFrame, delim=";")
    ########## Hydro data ##########
    xls_hydro = openxlsx(path*"hydro-v13-yw-2022_04_22.xlsx")
    Reservoir_ts = DataFrame(gettable(xls_hydro["reservoir"])...)
    RunOfRiver_ts = DataFrame(gettable(xls_hydro["run_of_river"])...)     
    ########## Hydrogen demand data ########## 
    Demand_hydrogen_df = DataFrame(gettable(
        openxlsx(path*"hydrogen-v02-yw-2022_04_25.xlsx")["demand"])...)
    ########## Heat demand data ##########
    min_load_ts_dict = 
        calculate_chp_min_load(copy(ES.countries), 2015, path, ES.T; offset=0.2)
    ########## Wind forecast data ########## 
    WindForecast_ts = CSV.read(joinpath(
        path, "scenarios_forecasterror-v04-yw-2022_05_10_flt.csv"), DataFrame, delim=";")
    dateformat = Dates.DateFormat("yyyy-mm-ddTHH:MM:SSZ")
    WindForecast_ts[!,:time] = Dates.DateTime.(WindForecast_ts[!,:time], dateformat)
    WindForecast_ts[!,:y] = Dates.year.(WindForecast_ts[!,:time])

    ES.forecast_error_scenarios = merge(
        [calculate_forecast_error_scenarios(
            WindForecast_ts,
            n,
            ES.T,
            ES.wind_uncertainty_scaling_factor) 
                for n in ES.N]...)

    RunOfRiver_ts, ES.hydro_sample_years = select_hydro_sample_years(
        RunOfRiver_ts, ES)

    ########## Caculate ramping rates conventionals ##########
    if any(occursin.("ramping_rates", readdir(path)))

        file_name = readdir(path)[occursin.("ramping_rates", readdir(path))][end]
        ramping_df = CSV.read(joinpath(path, file_name), DataFrame)
        ramping_rates = Dict()
        for c in ramping_df[!, :country]
            ramping_rates_c = Dict()
            for f in names(ramping_df)[names(ramping_df) .!= "country"]
                if ramping_df[ramping_df[!, :country] .== c, f][1] != 1.0
                    ramping_rates_c[f] = ramping_df[ramping_df[!, :country] .== c, f][1]
                end
            end
            ramping_rates[c] = ramping_rates_c
        end
        ES.ramping_rates = ramping_rates

    else
        ES.ramping_rates = 
            calculate_historical_technical_constraints(
                copy(ES.countries),
                ES.quantile_ramping,
                path*"historical_generation/")
    end

    Maxima_df = DataFrame(gettable(
        openxlsx(path*"maxima_2022_05_11.xlsx")["maxima"])...)

    ######################### Create Components #########################
    ############### DispatchableGenerators ###############
    Dispatchables_df = extract_tech_subset(Capacity_df, techs_dispatchable, ES)
    Conventionals = Dict()
    Reservoirs = Dict()
    for row in 1:nrow(Dispatchables_df)
        # Parameters
        id = Dispatchables_df[row, :n] * "_" * Dispatchables_df[row, :g]
        country, tech, g_max = Dispatchables_df[row, [:n, :g, :g_max]]
        fuel, eta, g_min_share, co2_fac, cost_vom = DataFrameRow(
            Parameters_i[Parameters_i[!, :i] .== tech,
                [:fuel, :avg_eta, :g_min_share, :co2_fac, :vom]], 1)
        g_min = repeat([g_max * g_min_share], length(ES.T))
        if occursin("chp", tech)
            g_min = g_max .* min_load_ts_dict[country]
        end
        if tech == "reservoir"
            g_min, g_max = 
                Maxima_df[
                    ((Maxima_df[!,:n] .== country) .&
                    (Maxima_df[!,:i] .== tech)), "ratio"][1] .* [g_min, g_max]
        end
        # Costs and prices
        cost_fuel = FuelPrice_df[
            ((FuelPrice_df[!,:fuel] .== fuel) .&
            (FuelPrice_df[!,:y] .== ES.year_model)), :p_fuel][1] / eta
        cost_CO2 = 
            co2_fac / eta *
                CO2Price_df[
                    ((CO2Price_df[!, :y] .== ES.year_model) .&
                     (CO2Price_df[!, :pathway] .== ES.co2_price_pathway)),
                        :p_co2][1]
        mc = cost_fuel + cost_CO2 + cost_vom
        bal_up, bal_down = mc .* [(1+ES.μ), (1-ES.μ)]

        r_up, r_down = nothing, nothing
        try
            if fuel in [
                "hardcoal", "biomass", "lignite", "nuclear", "gas", "peat"]
                r_up, r_down = repeat([g_max], 2) * ES.ramping_rates[country][fuel]
            elseif occursin("oil", fuel)
                r_up, r_down = repeat([g_max], 2) * ES.ramping_rates[country]["oil"]
            elseif fuel == "mix"
                fuels_availble = [f for f in keys(ES.ramping_rates[country])
                    if f in ["hardcoal", "lignite", "gas"]]
                # fuels_availble = keys(ES.ramping_rates[country])
                ramp_rate = minimum([ES.ramping_rates[country][f]
                    for f in fuels_availble])
                # Check that ramp does not conflict with chp minload
                r_up, r_down = repeat([g_max], 2) * 
                    ifelse(
                        maximum(abs.(diff(min_load_ts_dict[country]))) >= ramp_rate,
                        maximum(abs.(diff(min_load_ts_dict[country]))),
                        ramp_rate)
            elseif tech == "reservoir"
                r_up, r_down = repeat([g_max], 2) * ES.ramping_rates[country]["reservoir"]
            else
                r_up, r_down = repeat([g_max], 2)
            end
        catch KeyError
            println("Approximate ramping rate for fuel $fuel in country $country.")
            if occursin("oil", fuel)
                r_up, r_down = repeat([g_max], 2) * 
                mean([ES.ramping_rates[c]["oil"] for c in keys(ES.ramping_rates)
                    if "oil" in keys(ES.ramping_rates[c])])
            elseif ((fuel != "water") && (~occursin("oil", fuel)))
                r_up, r_down = repeat([g_max], 2) * 
                    mean([ES.ramping_rates[c][fuel] for c in keys(ES.ramping_rates)
                        if fuel in keys(ES.ramping_rates[c])])
            elseif tech == "reservoir"
                r_up, r_down = repeat([g_max], 2) * 
                    mean([ES.ramping_rates[c]["reservoir"] 
                        for c in keys(ES.ramping_rates)
                            if tech in keys(ES.ramping_rates[c])])
            else
                throw(KeyError)
            end
        end
        # Create components
        if ~(tech == "reservoir")
            Conventionals[id] = ConventionalGenerator(
                id, country, tech, g_min, g_max,
                bal_up, bal_down, r_up, r_down, mc
                )
        # Create extra set for reservoir plants to include generation limits
        else
            gen_limit = calculate_generation_limit(ES, Reservoir_ts, country)

            Reservoirs[id] = Reservoir(
                id, country, tech, g_min, g_max,
                bal_up, bal_down, r_up, r_down, mc, gen_limit
                )
            # Dispatchables[id] = Reservoirs[id]
        end

    end
    ES.Conventionals = Conventionals
    ES.Reservoirs = Reservoirs
    flt = [g == "reservoir" ? true : false for g in Dispatchables_df[!, :g]]
    ES.I = Dispatchables_df[.~flt, :n] .* "_" .* Dispatchables_df[.~flt, :g]
    ES.R = Dispatchables_df[flt, :n] .* "_" .* Dispatchables_df[flt, :g]

    # Create Superset of Dispatchables, Storages, and Electrolysers
    ES.U = vcat(ES.I, ES.R)
    ES.Units = merge(ES.Conventionals, ES.Reservoirs)


    ############### Storages ###############
    flt = [((s in techs_storage) && (y == ES.year_model) && (n in ES.countries)) ?
        true : false for (s,y,n) in eachrow(Storages_df[!, [:s, :y, :n]])]
    Storages_df = Storages_df[flt, :]
    Storages = Dict()
    for row in 1:nrow(Storages_df)
        # Parameters
        id = Storages_df[row, :n] * "_" * Storages_df[row, :s]
        country, tech, disc_max, ch_max, s_max = 
            Storages_df[row, [:n, :s, :disc_max, :ch_max, :e_max]]

        # Normalize charge and discharge rates for pumped hydro storage
        # Use mean of all others if no data exists
        if tech == "pumped_storage"
            try
                ch_max = Maxima_df[((Maxima_df[!,:n] .== country) .&
                    (Maxima_df[!,:i] .== "phs_charge")), "ratio"][1] * ch_max
                disc_max = Maxima_df[((Maxima_df[!,:n] .== country) .&
                    (Maxima_df[!,:i] .== "phs_discharge")), "ratio"][1] * disc_max
            catch BoundsError
                ch_max = mean(Maxima_df[
                    (Maxima_df[!,:i] .== "phs_charge"), "ratio"]) * ch_max
                disc_max = mean(Maxima_df[
                    (Maxima_df[!,:i] .== "phs_charge"), "ratio"]) * disc_max
            end
        end

        disc_eta, ch_eta, cost_vom, s_min_share, s_init_share =
            DataFrameRow(Parameters_s[Parameters_s[!, :s] .== tech,
                [:disc_eta, :ch_eta, :vom, :g_min_share, :e_init_share]], 1)
        ch_min, disc_min = s_min_share .* [ch_max, disc_max]
        s_init = s_init_share * s_max
        s_min = 0
        # Costs and prices
        mc = cost_vom

        r_up_L, r_down_L, r_up_G, r_down_G = repeat([nothing], 4)
        if tech == "pumped_storage"
            try
                r_up_L, r_down_L = repeat([ch_max], 2) * 
                    ES.ramping_rates[country]["phs_charge"]
                r_up_G, r_down_G = repeat([disc_max], 2) *
                    ES.ramping_rates[country]["phs_discharge"]
            catch KeyError
                println("Approximate ramping rate for $tech in country $country.")
                r_up_L, r_down_L = repeat([ch_max], 2) *
                    mean([ES.ramping_rates[c]["phs_charge"] for c in keys(ES.ramping_rates)
                        if "phs_charge" in keys(ES.ramping_rates[c])])
                r_up_G, r_down_G = repeat([disc_max], 2) *
                mean([ES.ramping_rates[c]["phs_discharge"] for c in keys(ES.ramping_rates)
                    if "phs_discharge" in keys(ES.ramping_rates[c])])
            end
        else
            r_up_L, r_down_L = repeat([ch_max], 2)
            r_up_G, r_down_G = repeat([disc_max], 2)
        end
            
        # Create components
        Storages[id] = Storage(
            id, country, tech, ch_min, ch_max, disc_min, disc_max,
            ch_eta, disc_eta, s_min, s_max, s_init,
            r_up_L, r_down_L, r_up_G, r_down_G, mc
            )
    end
    ES.Storages = Storages
    ES.SS = Storages_df[!, :n] .* "_" .* Storages_df[!, :s]

    ############### RenewableGenerators ###############
    Renewables_df = extract_tech_subset(Capacity_df, techs_res, ES)
    Renewables = Dict()
    for row in 1:nrow(Renewables_df)
        #Parameters
        id = Renewables_df[row, :n] * "_" * Renewables_df[row, :g] 
        country, tech, g_max = Renewables_df[row, [:n, :g, :g_max]]
        # Costs and prices
        cost_vom = Parameters_j[Parameters_j[!, :j] .== tech, :vom][1]
        mc = cost_vom

        ############### Overwriting Solar thermal ###############
        if tech == "solar_thermal"
            tech = "solarpv"
        end
        # Read in time series data
        if tech == "solarpv"
            point_forecast, scenario_output = 
                extract_ts_data_pv(ES, PV_ts, country, tech, g_max)
        elseif tech == "onshore_wind"
            point_forecast, scenario_output = 
                extract_ts_data_wind(ES, WindOnshore_ts, country, tech, g_max)
        elseif tech == "offshore_wind"
            point_forecast, scenario_output = 
                extract_ts_data_wind(ES, WindOffshore_ts, country, tech, g_max)
        elseif tech == "run_of_river"
            point_forecast, scenario_output = 
                extract_ts_data_ror(
                    ES,
                    RunOfRiver_ts,
                    country,
                    tech,
                    g_max)
        else
            throw(ArgumentError("RES technology is not in the predefined set."))
        end

        # Create components
        Renewables[id] = RenewableGenerator(
            id, country, tech, g_max, point_forecast,
            scenario_output, mc, nothing, nothing, nothing,
            nothing
            )
    end
    ES.Renewables = Renewables
    ES.J = Renewables_df[!, :n] .* "_" .* Renewables_df[!, :g]
    ES.Ω = isnothing(ES.years_scenarios) ? 
        ["s" * string(y-1979) for y in collect(1980:2019)] :
        ["s" * string(y) for y in 1:length(ES.years_scenarios)]

    ############### Electrolyzers ###############
    Electrolysers_df = extract_tech_subset(Capacity_df, techs_ptg, ES)
    Electrolysers = Dict()
    for row in 1:nrow(Electrolysers_df)
        # Parameters
        id = Electrolysers_df[row, :n] * "_" * Electrolysers_df[row, :g] 
        country, tech, l_max = Electrolysers_df[row, [:n, :g, :g_max]]
        eta, output_fuel, cost_vom = 
            DataFrameRow(Parameters_e[Parameters_e[!,:e] .== tech,
                [:eta, :fuel_output, :vom]], 1)

        # Costs and prices
        price_H2 = FuelPrice_df[
            ((FuelPrice_df[!,:fuel] .== output_fuel) .&
            (FuelPrice_df[!,:y] .== ES.year_model)), :p_fuel][1]
        mc = cost_vom

        p_B_up, p_B_down = mc .* [1-ES.μ, 1+ES.μ] 

        # Create components
        Electrolysers[id] = Electrolyser(
            id, country, tech, l_max,
            p_B_up, p_B_down, mc,
            price_H2, eta,
            nothing, nothing, Dict(), nothing
            )
    end
    ES.Electrolysers = Electrolysers
    ES.E = Electrolysers_df[!, :n] .* "_" .* Electrolysers_df[!, :g] 

    ############### Demands ###############
    Loads = Dict()
    for country in ES.countries
        if !(country in ES.countries)
            continue
        end
        demand_ts = []
        # Load time series data. Set demand to zero for nodes without demand (hubs)
        try
            demand_ts = Demand_EL_ts[Demand_EL_ts[!, :y] .== ES.year_model, country]
        catch e
            println("Set demand for $country to zero.")
            demand_ts = zeros(nrow(Demand_EL_ts))
        end
        # Value of lost load
        voll = 3000

        # Create components
        Loads[country] = Load(country, country, demand_ts, voll)
    end
    ES.Loads = Loads
    ES.D = ES.N # One demand per country

    ############### Lines ###############
    flt = [((n in ES.countries) && (m in ES.countries) && (y == ES.year_model)) ?
        true : false for (n,m,y) in eachrow(Interconnectors_df[!, [:n, :m, :y]])]
    Interconnectors_df = Interconnectors_df[flt, :]
    Interconnectors = Dict()
    for row in 1:nrow(Interconnectors_df)
        from, to, ntc = Interconnectors_df[row, [:n, :m, :ntc]]
        id = from * "_" * to
        # Create components
        Interconnectors[id] = Interconnector(id, from, to, ntc, nothing)
    end
    ES.Interconnectors = Interconnectors
    ES.Li = Interconnectors_df[!, :n] .* "_" .* Interconnectors_df[!, :m]

    # Ensure bottleneck free trade within bidding zone when "HBZ" is chosen
    if ES.bidding_zone_config == "HBZ"
        ES.Interconnectors["HUB1_DKW1"].ntc = 100000
        ES.Interconnectors["DKW1_HUB1"].ntc = 100000
        ES.Interconnectors["HUB3_DKE1"].ntc = 100000
        ES.Interconnectors["DKE1_HUB3"].ntc = 100000

        # for hub in ES.N[occursin.("HUB", ES.N)]
        #     for i in ES.Li[occursin.(hub, ES.Li)]
        #         ES.Interconnectors[i].ntc = 100000
        #     end
        # end
    end

    if ES.ntc_scaling_factor != 1.0
        for hub in ES.N[occursin.("HUB", ES.N)]
            for i in ES.Li[occursin.(hub, ES.Li)]
                ES.Interconnectors[i].ntc = 
                    ES.Interconnectors[i].ntc * ES.ntc_scaling_factor
            end
        end
    end



    ############### Hydrogen demand ###############
    flt = [((n in ES.countries) && (y == ES.year_model)) ?
        true : false for (n,y) in eachrow(Demand_hydrogen_df[!, [:n, :y]])]
    # Normalize the hydrogen with the amount of time steps considered
    ES.hydrogen_demand = (sum(Demand_hydrogen_df[flt, :h_tot])) * length(ES.T) / 8760;

end