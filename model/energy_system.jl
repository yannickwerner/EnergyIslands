mutable struct EnergySystem

    case_study # Name of case study containing the data files
    year_pointforecast # Year to use as point forecast for RES
    years_scenarios # List of weather years used as RES scenarios
    year_model # Model year (case study year)
    file_name # Name of excel file containing the data
    countries # List of countries to be modeled
    seed # Random seed for sampling of hydro timeseries data
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
    data_set_capacities # determines which data source to use
    ntc_scaling_factor # determines interconnector capacity between hubs and shore
    total_run_time # measures the total run time of a scenario including model building and solving
    scenario # name of the scenario

    # Dispatchables
    Conventionals
    Electrolysers
    Units
    Loads
    Interconnectors
    Renewables
    Demands
    Storages
    Reservoirs
    I
    D
    R
    T
    Li
    # G
    # D
    U
    J
    Ω
    E
    N
    # F
    SS

    ##### Solutions optimization model
    model_statistics
    G
    G_spill
    L
    L_shed
    F
    F_adj
    B_up
    B_down
    B_down_L
    B_up_L
    B_up_G
    B_down_G
    S
    S_RT
    excess
    shortage
    λ_DA
    λ_RT
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