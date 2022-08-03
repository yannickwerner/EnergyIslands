# Define Structs
mutable struct ConventionalGenerator
    id # id 
    node # bidding zone
    tech # generator technology
    g_min # minimum capacity
    g_max # maximum capacity
    p_B_up # bidding price (system) upward regulation
    p_B_down # bidding price (system) downward regulation
    r_up # ramp limit upward
    r_down # ramp limit downward
    mc # marginal cost of production
end


mutable struct Reservoir
    id # id 
    node # bidding zone
    tech # reservoir technology (hydro reservoir only)
    g_min # minimum capacity
    g_max # maximum capacity
    p_B_up # bidding price (system) upward regulation
    p_B_down # bidding price (system) downward regulation
    r_up # ramp limit upward
    r_down # ramp limit downward
    mc # marginal cost of production
    g_tot # maximum power production in given time period
end


mutable struct Storage
    id # id 
    node # bidding zone
    tech # storage technolgy (pumped-hydro or battery only)
    l_min # minimum capacity charging
    l_max # maximum capacity charging
    g_min # minimum capacity discharging
    g_max # maximum capacity discharging
    η_L # charging efficiency
    η_G # discharging efficiency
    s_min # minimum storage level
    s_max # maximum storage level
    s_init # initial storage level (default 50% of maximum)
    r_up_L # ramp limit upward charging mode
    r_down_L # ramp limit downward charging mode
    r_up_G # ramp limit upward discharging mode
    r_down_G # ramp limit downward discharging mode
    mc # marginal cost of production
end


mutable struct RenewableGenerator
    id # id 
    node # bidding zone
    tech # Renewable technology
    g_max # maximum capacity
    fc # probabilistic production forecast
    real # real-time realization power production
    mc # marginal cost of production
    expected_spill # ex-post calculation of expected spillage
    mean_expected_power_price # ex-post calculation of expected market value
    expected_profit # ex-post calculation of expected profit
    expected_production # ex-post calculation of total expected production
end

mutable struct Electrolyser
    id # id 
    node # bidding zone
    tech # electrolyser technology (onshore/offshore)
    l_max # maximum capacity (power)
    f_min # minimal capacity share (power, default 0)
    p_B_up # bidding price (system) upward regulation
    p_B_down # bidding price (system) downward regulation
    mc # marginal cost of production
    p_H2 # hydrogen price
    η # power-to-hydrogen efficiency
    profit_DA # ex-post calculation of day-ahead profits
    profit_BA # ex-post calculation of balancing market profits
    capacity_factors # ex-post calculation of capacity factors
    mean_expected_power_price # ex-post calculation of expected electricity costs
end

mutable struct Load
    id # id 
    node # bidding zone
    load # time-varying load
    voll # cost of load shedding (value of lost load)
end

mutable struct Interconnector
    id # id 
    from # starting bidding zone
    to # ending bidding zone
    ntc # maximum (net-transfer-)capacity
    expected_trade # ex-post calculation of expected trade
end