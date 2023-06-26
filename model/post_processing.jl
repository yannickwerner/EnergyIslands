function calculate_electrolyser_profit_DA(ES, e)
    """
    Calculates day-ahead profit of a given electrolyser
    based on the optimal solution of the energy system
    optimization model.

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    e : str
        Name of electrolyser

    Returns
    -------
    float
        Day-ahead profit of the given electrolyser
    """
    return sum(ES.L[e][t] * (
        ES.Electrolysers[e].η * ES.Electrolysers[e].p_H2-
        ES.Electrolysers[e].mc-
        ES.λ_DA[ES.Electrolysers[e].node][t]) for t in ES.T)
end


function calculate_electrolyser_profit_RT(ES, e)
    """
    Calculates balancing market profit of a given electrolyser
    based on the optimal solution of the energy system
    optimization model for every given scenario
    and on expectation.

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    e : str
        Name of electrolyser

    Returns
    -------
    float
        Probabilistic balancing market profit of the
        given electrolyser
    """
    profit_dict = Dict(ω => 
        sum(ES.B_down[e][ω][t] * (
                ES.Electrolysers[e].η * ES.Electrolysers[e].p_H2 -
                ES.Electrolysers[e].p_B_down -
                ES.λ_RT[ES.Electrolysers[e].node][ω][t]) -
            ES.B_up[e][ω][t] * (
                ES.Electrolysers[e].η * ES.Electrolysers[e].p_H2 -
                ES.Electrolysers[e].p_B_up -
                ES.λ_RT[ES.Electrolysers[e].node][ω][t])
            for t in ES.T) for ω in ES.Ω)
    profit_dict["expected"] = 1/length(ES.Ω) * sum(profit_dict[ω] for ω in ES.Ω)

    return profit_dict
end


function calculate_descriptive_statistics!(ES)
    """
    Calculates descriptive statistics for all electrolysers,
    renewable power plants, and interconnectors based on the
    optimal solutions from the energy system model
    including day-ahead profit, balancing market profit,
    capacity factors, and the mean expected electricity cost
    (electrolyser); expected spillage, expected power price,
    expected profit, and expected production.

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    Returns
    -------
    ES : EnergySystem
        Containing descriptive statistics
    """
    ############### Electrolyser ###############
    for e in ES.E
        ES.Electrolysers[e].profit_DA = 
            calculate_electrolyser_profit_DA(ES, e)
        ES.Electrolysers[e].profit_BA = 
            calculate_electrolyser_profit_RT(ES, e)
        ES.Electrolysers[e].capacity_factors["DA"] =
            mean(ES.L[e]) / ES.Electrolysers[e].l_max
        ES.Electrolysers[e].capacity_factors["BA_down"] = 
            mean(sum([ES.B_down[e][ω] for ω in ES.Ω]))/length(ES.Ω) / ES.Electrolysers[e].l_max
        ES.Electrolysers[e].capacity_factors["BA_up"] = 
            mean(sum([-ES.B_up[e][ω] for ω in ES.Ω]))/length(ES.Ω) / ES.Electrolysers[e].l_max
        ES.Electrolysers[e].mean_expected_power_price = 
            sum(sum(1/length(ES.Ω) *
                (ES.L[e]/(1/length(ES.Ω))/length(ES.Ω) .* ES.λ_DA[ES.Electrolysers[e].node] + 
                 ES.B_down[e][ω]-ES.B_up[e][ω] .* ES.λ_RT[ES.Electrolysers[e].node][ω]
                ) for ω in ES.Ω)) /
            1/sum((sum(1/length(ES.Ω) *
            (ES.L[e]/(1/length(ES.Ω))/length(ES.Ω) + 
             ES.B_down[e][ω] - ES.B_up[e][ω]
            ) for ω in ES.Ω)))
    end

    ############### Renewables ###############
    for j in ES.J
        ES.Renewables[j].expected_spill =
            vec(mean(Array(DataFrame(ES.G_spill[j])), dims=2))
        ES.Renewables[j].mean_expected_power_price =
            sum(sum(1/length(ES.Ω) *
                (ES.G[j]/(1/length(ES.Ω))/length(ES.Ω) .* ES.λ_DA[ES.Renewables[j].node] + 
                (ES.Renewables[j].real[ω]-ES.G[j]-ES.G_spill[j][ω]) .* ES.λ_RT[ES.Renewables[j].node][ω]
                ) for ω in ES.Ω)) /
            1/sum((sum(1/length(ES.Ω) *
                (ES.G[j]/(1/length(ES.Ω))/length(ES.Ω) + 
                (ES.Renewables[j].real[ω]-ES.G[j]-ES.G_spill[j][ω])
                ) for ω in ES.Ω)))
        ES.Renewables[j].expected_profit =
            sum(ES.G[j] .* ES.λ_DA[ES.Renewables[j].node] +
            1/length(ES.Ω) * 
                sum((ES.Renewables[j].real[ω]-ES.G[j]-ES.G_spill[j][ω]) .*
                ES.λ_RT[ES.Renewables[j].node][ω]
                for ω in ES.Ω))
        ES.Renewables[j].expected_production = 
            mean([ES.Renewables[j].real[ω] for ω in ES.Ω])
    end

    ############### Interconnectors ###############
    for l in ES.Li
        ES.Interconnectors[l].expected_trade = 
            ES.F[l] .+ vec(mean(Array(DataFrame(ES.F_adj[l])), dims=2))
    end

end


function calculate_non_local_renewable_power_consumption_electrolyers(ES, e, ω)
    """
    Calculates and returns the non-local renewable power consumption of a given electrolyzer
    in a given scenario.

    Parameters
    ----------    
    ES : EnergySystem
        Contains information about energy system modeled

    e : str
        Name of electrolyser

    ω : str
        Name of scenario

    Returns
    -------
    non_local_renewable_power_consumption : Vector{Float64}
        Non-local renewable power consumption of electrolyzer e in scenario ω
    """
    local_renewables = [j for j in ES.J if ES.Renewables[j].node ==  ES.Electrolysers[e].node]
    local_production_real = reduce(+, [
        ES.Renewables[j].real[ω] - ES.G_spill[j][ω] for j in local_renewables])
    non_local_renewable_power_consumption = 
        ES.L[e] + ES.B_down[e][ω] - ES.B_up[e][ω]  - local_production_real
    non_local_renewable_power_consumption[non_local_renewable_power_consumption .<= 0] .= 0

    return non_local_renewable_power_consumption
end


function extract_optimal_values(scenarios, countries, ES)
    """
    Extracts specific optimal solutions for certain variables
    from the model in order to save them later on.

    Parameters
    ----------
    scenarios : Array
        Contains scenarios for which to retrieve optimal
        solutions from the model.

    countries : Array
        Contains countries for which to retrieve optimal
        soltuions from the model
    
    ES : EnergySystem
        Contains information about energy system modeled

    Returns
    -------
    df : DataFrame
        Contains optimal solutions from the model
    """
    df = DataFrame()
    variable_list_DA = [:F, :L, :G, :λ_DA]
    for var in variable_list_DA
        df_sub = DataFrame(getproperty(ES, var))
        rename!(df_sub, string(var)*"_".*names(df_sub))
        df = hcat(df, df_sub)
    end

    variable_list_RT = [
        :F_adj, :G_spill, :L_shed, :B_down, :B_up, :B_down_L,
        :B_up_L, :B_down_G, :B_up_G, :λ_RT]
    for var in variable_list_RT
        data = getproperty(ES, var)
        var_keys = collect(keys(data))
        for key in var_keys
            df_sub = DataFrame(data[key])
            var_expected = reduce(+, eachcol(df_sub)) ./ ncol(df_sub)
            select!(df_sub, scenarios)
            rename!(df_sub, string(var)*"_"*key*"_".*names(df_sub))
            df = hcat(df, df_sub)
            df = hcat(df, DataFrame(
                Symbol(string(var)*"_"*key*"_avg") => var_expected))
        end
    end

    cols_select = Bool.(zeros(length(names(df))))
    for c in countries
        cols_select = cols_select .| occursin.(c, names(df))
    end
    select!(df, cols_select)

    col_names = names(df)
    for col in 1:length(col_names)
        if occursin.("λ", col_names[col])
            col_names[col] = replace(col_names[col], "λ" => "price")
        end
    end
    rename!(df, col_names)

    for j in [j for j in ES.J if 
        ES.Renewables[j].node in countries]
        for ω in scenarios
            df[!, "Deviation_$(j)_$(ω)"] =
                ES.Renewables[j].real[ω] .- ES.G[j]
        end
    end

    # Electrolyser electricity imports
    for e in ES.E
        nl_power_cons_e_sc = reduce(hcat,
            calculate_non_local_renewable_power_consumption_electrolyers.(
                Ref(ES), Ref(e), ES.Ω))
        nl_power_cons_e = reduce(vcat, mean(nl_power_cons_e_sc, dims=2))
        df[!, "nl_res_power_cons_avg_$e"] = nl_power_cons_e
    end


    # System data
    for ω in scenarios
        df_system = DataFrame(
            "system_spillage_RES_$ω" =>
                sum([ES.G_spill[j][ω][ES.T] for j in ES.J]),
            "system_deviation_RES_$ω" =>
                sum([ES.Renewables[j].real[ω][ES.T].-ES.G[j][ES.T] for j in ES.J]),
            "system_curtailment_load_$ω" => 
                sum([ES.L_shed[n][ω][ES.T] for n in ES.N]),
            "system_balancing_electrolyser_$ω" =>
                sum([ES.B_down[e][ω][ES.T] -
                    ES.B_up[e][ω][ES.T] for e in ES.E]),
            "system_balancing_storage_$ω" =>
                sum([ES.B_down_L[s][ω][ES.T] + ES.B_down_G[s][ω][ES.T] - (
                    ES.B_up_L[s][ω][ES.T] + ES.B_up_G[s][ω][ES.T])
                    for s in ES.SS]),
            "system_balancing_reservoir_$ω" =>
                sum([ES.B_down[r][ω][ES.T] -
                    ES.B_up[r][ω][ES.T] for r in ES.R]),
            "system_balancing_conventionals_$ω" =>
                sum([ES.B_down[i][ω][ES.T] - 
                    ES.B_up[i][ω][ES.T] for i in ES.I])
        )
        
        df = hcat(df, df_system)
    end

    ##### Storage simulatenous balancing #####
    activations = zeros(length(ES.T))
    magnitude = zeros(length(ES.T))
    for ss in ES.SS
        for ω in ES.Ω
            up_activation =
                (ES.L[ss] .> 0.0) .&
                (ES.B_up_L[ss][ω] .> 0.0) .&
                (ES.B_up_G[ss][ω] .> 0.0)
            magnitude += 
                up_activation .* ES.B_up_G[ss][ω]

            down_activation =
                (ES.G[ss] .> 0.0) .&
                (ES.B_down_G[ss][ω] .> 0.0) .&
                (ES.B_down_L[ss][ω] .> 0.0)
            magnitude += 
                down_activation .* ES.B_down_L[ss][ω]
            
            activations += 
                float.(up_activation) +
                float.(down_activation)
        end
    end

    df[!, "stor_double_bal_activations"] = activations
    df[!, "stor_double_bal_magnitude"] = magnitude

    df = round.(df, digits=2)
    df = Float32.(df)

    return df

end


function calculate_objective_value_scenario(ES)
    """
    Calculates total system costs (objective value) for each scenario.

    Parameters
    ----------
     ES : EnergySystem
        Contains information about energy system modeled

    Returns
    -------
    obj_val : Dict
        Contains total system costs (objective value) for 
        each scenario
    """
    obj_val = Dict()
    for ω in ES.Ω
        if ES.hydrogen_production_driver == "price"
            obj_val[ω] = sum(
                sum(ES.Conventionals[i].mc*ES.G[i][t] for i in ES.I) +
                sum((-ES.Electrolysers[e].p_H2*ES.Electrolysers[e].η+
                    ES.Electrolysers[e].mc)*ES.L[e][t] for e in ES.E) +
                sum(ES.Conventionals[i].p_B_up*ES.B_up[i][ω][t] - 
                    ES.Conventionals[i].p_B_down*ES.B_down[i][ω][t] for i in ES.I) +
                sum((ES.Electrolysers[e].p_H2*ES.Electrolysers[e].η-
                    ES.Electrolysers[e].p_B_up)*ES.B_up[e][ω][t] +
                    (-ES.Electrolysers[e].p_H2*ES.Electrolysers[e].η+
                        ES.Electrolysers[e].p_B_down)*ES.B_down[e][ω][t] for e in ES.E) +
                sum(ES.Loads[l].voll*ES.L_shed[l][ω][t] for l in ES.D) 
                for t in ES.T)
                    
        elseif ES.hydrogen_production_driver == "demand"
            obj_val[ω] = sum(
                sum(ES.Conventionals[i].mc*ES.G[i][t] for i in ES.I) +
                sum(Electrolysers[e].mc*ES.L[e][t] for e in ES.E) +
                sum(ES.Conventionals[i].p_B_up*ES.B_up[i][ω][t] - 
                    ES.Conventionals[i].p_B_down*ES.B_down[i][ω][t] for i in ES.I) +
                sum(-ES.Electrolysers[e].p_B_up*ES.B_up[e][ω][t] +
                    ES.Electrolysers[e].p_B_down*ES.B_down[e][ω][t] for e in ES.E) +
                sum(ES.Loads[l].voll*ES.L_shed[l][ω][t] for l in ES.D) 
                for t in ES.T)
        end
    end
    return obj_val
end


function write_descriptive_statistics(ES::EnergySystem, weeks)
    """
    Retrieves and calculates descriptive statistics for the optimization
    models.

    Parameters
    ----------    
    ES : EnergySystem
        Contains information about energy system modeled

    weeks : Array
        Contains tuples of weeks to run the model for different
        time horizons

    Returns
    -------
    electrolyser_df : DataFrame
        Contains descriptive statistics of electrolysers
    renewables_df : DataFrame
        Contains descriptive statistics of renewable power plants
    node_df : DataFrame
        Contains descriptive statistics of bidding zones
    model_statistics_df : DataFrame
        Contains descriptive statistics of the model
    objective_value_df : DataFrame
        Contains objective values for each model time horion

    """

    hour_start = (weeks[1]-1)*168+1
    hour_end = weeks[end]*168

    electrolyser_df = DataFrame()
    renewables_df = DataFrame()
    node_df = DataFrame()

    for e in ES.E

        nl_power_cons_e_sc = reduce(hcat,
            calculate_non_local_renewable_power_consumption_electrolyers.(
                Ref(ES), Ref(e), ES.Ω))
        nl_power_cons_e = reduce(vcat, mean(nl_power_cons_e_sc, dims=2))

        df = DataFrame(
            :scenario => ES.scenario,
            :bz_config => ES.bidding_zone_config,
            :y => ES.year_model,
            :data_set => ES.data_set_capacities,
            :ntc_scaling_factor => ES.ntc_scaling_factor,
            :e => e,
            :n => ES.Electrolysers[e].node,
            :weeks => string(weeks),
            :hour_start => hour_start,
            :hour_end => hour_end,
            :no_scenarios => length(ES.Ω),
            :mean_expected_power_price =>
                ES.Electrolysers[e].mean_expected_power_price,
            :profit_DA => ES.Electrolysers[e].profit_DA,
            :profit_BA_expected => ES.Electrolysers[e].profit_BA["expected"],
            :profit_total => 
                ES.Electrolysers[e].profit_DA+
                ES.Electrolysers[e].profit_BA["expected"],
            :non_local_renewable_electricity_consumption =>
                sum(nl_power_cons_e)/(ES.Electrolysers[e].l_max * length(ES.T))
        )

        df = hcat(df, DataFrame(ES.Electrolysers[e].capacity_factors))

        electrolyser_df = reduce(vcat, [electrolyser_df, df])
    end 

    for j in ES.J
        df = DataFrame(
            :scenario => ES.scenario,
            :bz_config => ES.bidding_zone_config,
            :y => ES.year_model,
            :data_set => ES.data_set_capacities,
            :ntc_scaling_factor => ES.ntc_scaling_factor,
            :j => j,
            :n => ES.Renewables[j].node,
            :weeks => string(weeks),
            :hour_start => hour_start,
            :hour_end => hour_end,
            :no_scenarios => length(ES.Ω),
            :mean_expected_power_price =>
                ES.Renewables[j].mean_expected_power_price,
            :expected_spill => sum(ES.Renewables[j].expected_spill),
            :expected_profit => ES.Renewables[j].expected_profit
        )
        renewables_df = reduce(vcat, [renewables_df, df])
    end

    for n in ES.N
        df = DataFrame(
            :scenario => ES.scenario,
            :bz_config => ES.bidding_zone_config,
            :weeks => string(weeks),
            :hour_start => hour_start,
            :hour_end => hour_end,
            :y => ES.year_model,
            :ntc_scaling_factor => ES.ntc_scaling_factor,
            :data_set => ES.data_set_capacities,
            :no_scenarios => length(ES.Ω),
            :n => n,
            :mean_power_price_DA => mean(ES.λ_DA[n])
        )
        node_df = reduce(vcat, [node_df, df])
    end

    model_statistics_df = DataFrame(ES.model_statistics)
    model_statistics_df[!, :total_run_time] = [ES.total_run_time]
    model_statistics_df[!, :weeks] = [string(weeks)]
    model_statistics_df[!, :hour_start] = [hour_start]
    model_statistics_df[!, :hour_end] = [hour_end]
    model_statistics_df[!, :ntc_scaling_factor] =
        [ES.ntc_scaling_factor]

    objective_value_df = DataFrame(
        :scenario => ES.scenario,
        :bz_config => ES.bidding_zone_config,
        :y => ES.year_model,
        :data_set => ES.data_set_capacities,
        :ntc_scaling_factor => ES.ntc_scaling_factor,
        :weeks => string(weeks),
        :hour_start => hour_start,
        :hour_end => hour_end,
        :no_scenarios => length(ES.Ω))
    objective_value_df = hcat(
        objective_value_df,
        DataFrame(calculate_objective_value_scenario(ES))
    )

    return electrolyser_df, renewables_df, node_df, model_statistics_df, objective_value_df

end