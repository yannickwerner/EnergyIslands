mutable struct EnergySystem

    case_study # Name of case study containing the data files
    year_pointforecast # Year to use as point forecast for RES
    years_scenarios # List of weather years used as RES scenarios
    year_model # Model year (case study year)
    file_name # Name of excel file containing the data
    countries # List of countries to be modeled
    # seed # Random seed for sampling of hydro timeseries data
    hydro_sample_years # collects the randomly sampled years for hydro
    μ # cost margin for bid prices in balancing market
    hydrogen_demand # total hydrogen demand per bidding zone
    hydrogen_autarky_rate # share of hydrogen production from modeled countries
    forecast_error_scenarios # scenarios of forecast errors for RES
    bidding_zone_config # hub bidding zones; "OBZ" or "HBZ"
    model # container for the JuMP model
    quantile_ramping # quantile to use for calculating historical ramping rates
    ramping_rates # fuel-specific ramping rates based on historical data
    hydrogen_production_driver # includes hydrogen prices in objective function or hydrogen demand constraint
    uncertainty_horizon # wind power production forecast horizon
    wind_uncertainty_scaling_factor # factor used to scale wind uncertainty
    co2_price_pathway # determines CO2 prices
    h2_price_pathway # determines H2 prices
    data_set_capacities # determines which data source to use
    ntc_scaling_factor # determines interconnector capacity between hubs and shore
    total_run_time # measures the total run time of a scenario including model building and solving
    scenario # name of the scenario

    ##### Component data #####
    Conventionals # Conventional generator data
    Electrolysers # Electrolyzer data
    Units # Union of conventionals and electrolyzers
    Renewables # Renewable generator data
    Loads # Load data
    Interconnectors # Interconnector data
    #Demands
    Reservoirs # Hydro reservoir data
    Storages # Storage data

    ##### Sets #####
    I # Set of conventionals
    E # Set of electrolyzers
    U # Super-set of conventionals and reservoirs
    J # Set of renewables
    D # Set of loads
    Li # Set of interconnectors
    R # Set of reservoirs
    SS # Set of storages

    ## G
    ## D
    ## F

    T # Set of time steps
    Ω # Set of scenarios
    N # Set of nodes

    ##### Solutions optimization model #####
    model_statistics # Solution and model statistics
    G # Day-ahead schedules of dispatchable generators (I,E,J,R,SS) 
    G_spill # Real-time spillage of renewables (R) 
    L # Day-ahead schedule of dispatchable loads (D, SS)
    L_shed # Real-time load shedding (D)
    F # Day-ahead interconnector flows (Li)
    F_adj # Real-time adjustment of interconnector flows (Li)
    B_up # Real-time upward balancing of dispatchable generators (I,E,R)
    B_down # Real-time downwrad balancing of dispatchable generators (I,E,R)
    B_down_L # Real-time downward balancing of storages (SS)
    B_up_L # Real-time upward balancing of storages (SS)
    B_up_G # Real-time upward balancing of storages (SS)
    B_down_G # Real-time downward balancing of storages (SS)
    S # Day-ahead energy storage level (SS)
    S_RT # Real-time energy storage level (SS)
    excess # Real-time excess sink to ensure feasibility
    shortage # Real-time shortage source to ensure feasibility
    λ_DA # Day-ahead nodal prices
    λ_RT # Real-time nodal prices
    ρ_H2 # endogenous hydrogen price
    γ_H2 # endogenous hydrogen demand

    function EnergySystem(kwargs...)
        ES = new()
        for (key, value) in kwargs
            setfield!(ES, key, value)
        end
        return ES
    end
end