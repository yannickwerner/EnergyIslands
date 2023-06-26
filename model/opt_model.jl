function run_opt_model!(ES; print_model=false)
    """
    Builds and solves JuMP optimization model based on
    the given energy system. Obtains and returns the
    optimal solutions from the model, and writes those
    on the EnergySystem object passed.

    Parameters
    ----------
    ES : EnergySystem
        Contains information about energy system modeled

    print_model : str
        Option to choose if model is printed (true/false).
        Default false

    Returns
    -------
    ES : EnergySystem
        Contains information about energy system modeled
        and holds optimal solutions of the optimization
        problem.
    """

    # Extract sets and parameters from EnergySystem object
    Electrolysers = ES.Electrolysers
    Units = ES.Units
    Loads = ES.Loads
    Interconnectors = ES.Interconnectors
    Renewables = ES.Renewables
    Storages = ES.Storages
    Reservoirs = ES.Reservoirs
    Conventionals = ES.Conventionals
    T = ES.T
    SS = ES.SS
    Li = ES.Li
    U = ES.U
    J = ES.J
    Ω = ES.Ω
    E = ES.E
    N = ES.N
    R = ES.R
    D = ES.D
    I = ES.I

    T_red = T[2:end]

    B = vcat(ES.I, ES.R, ES.SS, ES.J)
    P = vcat(ES.I, ES.R, ES.E)
    A = vcat(ES.E, ES.SS)

    m = Model(Gurobi.Optimizer)
    # set_optimizer_attribute(m, "DualReductions", 0)
    # set_optimizer_attribute(m, "Threads", 24)
    # set_optimizer_attribute(m, "Method", 1) # uses dual simplex

    @variables m begin
        G[g in B, t in T] >= 0
        L[a in A, t in T] >= 0
        B_up[p in P, t in T, ω in Ω] >= 0
        B_down[p in P, t in T, ω in Ω] >= 0
        B_aux[u in U, t in T, ω in Ω]
        0 <= G_spill[j in J, t in T, ω in Ω] <= Renewables[j].real[ω][t]
        0 <= L_shed[d in D, t in T, ω in Ω] <= Loads[d].load[t]
        0 <= F[l in Li, t in T] <= Interconnectors[l].ntc
        0 <= F_adj[l in Li, t in T, ω in Ω] <= Interconnectors[l].ntc
        Storages[s].s_min <= S_RT[s in SS, t in T, ω in Ω] <= Storages[s].s_max
        Storages[s].s_min <= S[s in SS, t in T] <= Storages[s].s_max
        B_down_L[s in SS, t in T, ω in Ω] >= 0
        B_up_L[s in SS, t in T, ω in Ω] >= 0
        B_up_G[s in SS, t in T, ω in Ω] >= 0
        B_down_G[s in SS, t in T, ω in Ω] >= 0
        B_G_aux[s in SS, t in T, ω in Ω]
        B_L_aux[s in SS, t in T, ω in Ω]
    end
    
    
    if ES.hydrogen_production_driver == "price"
        @objective(m, Min,
            sum(
                sum(Conventionals[i].mc*G[i,t] for i in I) +
                sum(G[s,t]*0.001 for s in SS) +
                sum(F[l,t]*0.001 for l in Li) +
                sum((-Electrolysers[e].p_H2*Electrolysers[e].η+Electrolysers[e].mc)*L[e,t] for e in E)
                for t in T) +
            sum(1/length(Ω)*(sum(
                sum(Conventionals[i].p_B_up*B_up[i,t,ω] - 
                    Conventionals[i].p_B_down*B_down[i,t,ω] for i in I) +
                sum((Electrolysers[e].p_H2*Electrolysers[e].η-Electrolysers[e].p_B_up)*B_up[e,t,ω] +
                    (-Electrolysers[e].p_H2*Electrolysers[e].η+Electrolysers[e].p_B_down)*B_down[e,t,ω] for e in E) +
                sum(B_up_L[s,t,ω]*0.001 + B_up_G[s,t,ω]*0.001 for s in SS) +
                sum(F_adj[l,t,ω]*0.001 for l in Li) +
                sum(Loads[l].voll*L_shed[l,t,ω] for l in D) 
                for ω in Ω))
                for t in T))
    elseif ES.hydrogen_production_driver == "demand"
        @objective(m, Min,
        sum(
            sum(Conventionals[i].mc*G[i,t] for i in I) +
            sum(G[s,t]*0.001 for s in SS) +
            sum(F[l,t]*0.001 for l in Li) +
            sum(Electrolysers[e].mc*L[e,t] for e in E)
            for t in T) +
        sum(1/length(Ω)*(sum(
            sum(Conventionals[i].p_B_up*B_up[i,t,ω] - 
                Conventionals[i].p_B_down*B_down[i,t,ω] for i in I) +
            sum(-Electrolysers[e].p_B_up*B_up[e,t,ω] +
                Electrolysers[e].p_B_down*B_down[e,t,ω] for e in E) +
            sum(B_up_L[s,t,ω]*0.001 + B_up_G[s,t,ω]*0.001 for s in SS) +
            sum(F_adj[l,t,ω]*0.001 for l in Li) +
            sum(Loads[l].voll*L_shed[l,t,ω] for l in D) 
            for ω in Ω))
            for t in T))
    else
        throw(ArgumentError("Hydrogen production driver is set
            to $(ES.hydrogen_production_driver), but must be either
            'price' or 'demand'."))
    end

    @constraints m begin
        ########## Supply-demand balances ##########
        # Nodebalance day ahead
        NodeBalanceDA[n in N, t in T],
            sum(G[u,t] for u in U if Units[u].node == n) +
            sum(G[s,t] for s in SS if Storages[s].node == n) +
            sum(G[j,t] for j in J if Renewables[j].node == n) -
            sum(Loads[l].load[t] for l in D if Loads[l].node == n) -
            sum(L[e,t] for e in E if Electrolysers[e].node == n) -
            sum(L[s,t] for s in SS if Storages[s].node == n) -
            sum(F[f,t] for f in Li if Interconnectors[f].from == n) +
            sum(F[f,t] for f in Li if Interconnectors[f].to == n) == 0

        # Nodebalance real time
        NodeBalanceRT[n in N, t in T, ω in Ω],
            sum(Renewables[j].real[ω][t] - G[j,t] - G_spill[j,t,ω] for j in J if Renewables[j].node == n) +
            sum(B_up[u,t,ω] - B_down[u,t,ω] for u in U if Units[u].node == n) +
            sum(B_up[e,t,ω] - B_down[e,t,ω] for e in E if Electrolysers[e].node == n) +
            sum(B_up_G[s,t,ω] + B_up_L[s,t,ω] - B_down_G[s,t,ω] - B_down_L[s,t,ω] for s in SS if Storages[s].node == n) +
            sum(L_shed[d,t,ω] for d in D if Loads[d].node == n) -
            sum(F_adj[f,t,ω] for f in Li if Interconnectors[f].from == n) +
            sum(F_adj[f,t,ω] for f in Li if Interconnectors[f].to == n) == 0

        ########## Capacity constraints for conventional reservoir units ##########
        # Renewable constraints
        RenewablesLimitUp[j in J, t in T],
            G[j,t] <= Renewables[j].g_max
        # Generation limit dispachable gnerators
        GenerationLimitUpDispatchables[u in U, t in T, ω in Ω],
            G[u,t] + B_up[u,t,ω] <= Units[u].g_max
        GenerationLimitDownDispatchables[u in U, t in T, ω in Ω],
            G[u,t] - B_down[u,t,ω] >= Units[u].g_min[t]
        # Reservoir
        ReservoirTotalGenerationLimit[r in R, ω in Ω],
            sum(G[r,t]+B_up[r,t,ω]-B_down[r,t,ω] for t in T) <= Reservoirs[r].g_tot[ω]

        # Ramping constraints
        RampingLimitUpDispatchablesDA[u in U, t in T_red],
            G[u,t] - G[u,t-1] <= Units[u].r_up
        RampingLimitDownDispatchablesDA[u in U, t in T_red],
            G[u,t] - G[u,t-1] >= -Units[u].r_down
        RampingAuxiliaryVariableDefinitionDispatchablesRT[u in U, t in T, ω in Ω],
            B_aux[u,t,ω] == B_up[u,t,ω] - B_down[u,t,ω]
        RampingLimitUpDispatchablesRT[u in U, t in T_red, ω in Ω],
            G[u,t] - G[u,t-1] + B_aux[u,t,ω] -
                B_aux[u,t-1,ω] <= Units[u].r_up
        RampingLimitDownDispatchablesRT[u in U, t in T_red, ω in Ω],
            G[u,t] - G[u,t-1] + B_aux[u,t,ω] - 
                B_aux[u,t-1,ω] >= -Units[u].r_down

        # Interconnector flow limits
        FlowLimitNTCRT[l in Li, t in T, ω in Ω],
            F[l,t] + F_adj[l,t,ω] <= Interconnectors[l].ntc

        ############### Electrolyser Constraints ###############
        DemandLimitUpElectrolyser[e in E, t in T, ω in Ω],
            L[e,t] + B_down[e,t,ω] <= Electrolysers[e].l_max
        DemandLimitDownElectrolyzer[e in E, t in T, ω in Ω],
            L[e,t] - B_up[e,t,ω] >= Electrolysers[e].f_min * Electrolysers[e].l_max

        ############### Storage Constraints ###############
        ##### Day-ahead constraints #####
        ChargeEnergyLimitUp[s in SS, t in T, ω in Ω],
            L[s,t] + B_down_L[s,t,ω] <= Storages[s].l_max
        ChargeEnergyLimitDown[s in SS, t in T, ω in Ω],
            L[s,t] - B_up_L[s,t,ω] >= Storages[s].l_min
        DischargeEnergyLimitUp[s in SS, t in T, ω in Ω],
            G[s,t] + B_up_G[s,t,ω] <= Storages[s].g_max
        DischargeEnergyLimitDown[s in SS, t in T, ω in Ω],
            G[s,t] - B_down_G[s,t,ω] >= Storages[s].g_min

        ##### Ramping constraints #####
        RampingLimitUpStoragesDischargeDA[s in SS, t in T_red],
            G[s,t] - G[s,t-1] <= Storages[s].r_up_G
        RampingLimitDownStoragesDischargeDA[s in SS, t in T_red],
            G[s,t] - G[s,t-1] >= -Storages[s].r_down_G
        RampingLimitUpStoragesChargeDA[s in SS, t in T_red],
            L[s,t] - L[s,t-1] <= Storages[s].r_up_L
        RampingLimitDownStoragesChargeDA[s in SS, t in T_red],
            L[s,t] - L[s,t-1] >= -Storages[s].r_down_L

        RampingAuxiliaryVariableLDefinitionStoragesRT[s in SS, t in T, ω in Ω],
            B_L_aux[s,t,ω] == B_up_L[s,t,ω] - B_down_L[s,t,ω]
        RampingAuxiliaryVariableGDefinitionStoragesRT[s in SS, t in T, ω in Ω],
            B_G_aux[s,t,ω] == B_up_G[s,t,ω] - B_down_G[s,t,ω]

        RampingLimitUpStoragesDishargeRT[s in SS, t in T_red, ω in Ω],
            G[s,t] - G[s,t-1] + B_G_aux[s,t,ω] - B_G_aux[s,t-1,ω] <= Storages[s].r_up_G
        RampingLimitDownStoragesDishargeRT[s in SS, t in T_red, ω in Ω],
            G[s,t] - G[s,t-1] + B_G_aux[s,t,ω] - B_G_aux[s,t-1,ω] >= -Storages[s].r_down_G
        RampingLimitUpStoragesChargeRT[s in SS, t in T_red, ω in Ω],
            L[s,t] - L[s,t-1] + B_L_aux[s,t,ω] - B_L_aux[s,t-1,ω] <= Storages[s].r_up_L
        RampingLimitDownStoragesChargeRT[s in SS, t in T_red, ω in Ω],
            L[s,t] - L[s,t-1] + B_L_aux[s,t,ω] - B_L_aux[s,t-1,ω] >= -Storages[s].r_down_L

        ##### Storage level constraints #####
        StorageLevelDA[s in SS, t in T_red],
            S[s,t] == S[s,t-1] +
                Storages[s].η_L * L[s,t] - 1/Storages[s].η_G * G[s,t]
        StorageLevelInitDA[s in SS, t = [T[1]]],
            S[s,t] == Storages[s].s_init +
                Storages[s].η_L * L[s,t] - 1/Storages[s].η_G * G[s,t]
        StorageLevelEndDA[s in SS, t = [T[end]]],
            S[s,t] >= Storages[s].s_init

        StorageLevelRT[s in SS, t in T_red, ω in Ω],
            S_RT[s,t,ω] == S_RT[s,t-1,ω] +
                Storages[s].η_L * (L[s,t] + B_down_L[s,t,ω] - B_up_L[s,t,ω]) -  
                1/Storages[s].η_G * (G[s,t] - B_down_G[s,t,ω] + B_up_G[s,t,ω])
        StorageLevelInitRT[s in SS, t = [T[1]], ω in Ω],
            S_RT[s,t,ω] == Storages[s].s_init +
                Storages[s].η_L * (L[s,t] + B_down_L[s,t,ω] - B_up_L[s,t,ω]) -  
                1/Storages[s].η_G * (G[s,t] - B_down_G[s,t,ω] + B_up_G[s,t,ω])
        StorageLevelEndRT[s in SS, t = [T[end]], ω in Ω],
            S_RT[s,t,ω] >= Storages[s].s_init
    end
    ############### Hydrogen demand constraint ###############
    if ES.hydrogen_production_driver == "demand"
        @constraint(m, HydrogenDemand,
            sum(Electrolysers[e].η*(L[e,t] + 
                1/length(Ω)*sum(B_down[e,t,ω]-B_up[e,t,ω] for ω in Ω))
                for e in E, t in T) == 
                    ES.hydrogen_demand*ES.hydrogen_autarky_rate)
    end

    if print_model == true
        print(m)
    end

    @time optimize!(m)
    status = termination_status(m)
    println(status)
    println(raw_status(m))

    ES.model_statistics = Dict(
        :scenario => ES.scenario,
        :bidding_zone_config => ES.bidding_zone_config,
        :model_year => ES.year_model,
        :data_set_capacities => ES.data_set_capacities,
        :co2_price_pathway => ES.co2_price_pathway,
        :hours => length(ES.T),
        :ntc_scaling_factor => ES.ntc_scaling_factor,
        :pf_res_year => ES.year_pointforecast,
        :no_scenarios => length(ES.years_scenarios), 
        :objective_value => objective_value(m),
        :solution_time => solve_time(m),
        :number_constraints => sum(
            num_constraints(m, F, S) for (F, S) in
            list_of_constraint_types(m)),
        :number_variables => length(all_variables(m))
    )

    ES.G = Dict(b => [value(G[b,t]) for t in T] for b in B)
    ES.G_spill = Dict(j => Dict(ω => [value(G_spill[j,t,ω]) for t in T] for ω in Ω) for j in J)

    ES.L_shed = Dict(d => Dict(ω => [value(L_shed[d,t,ω]) for t in T] for ω in Ω) for d in D)

    ES.L = Dict(a => [value(L[a,t]) for t in T] for a in A)
    ES.S_RT = Dict(s => Dict(ω => [value(S_RT[s,t,ω]) for t in T] for ω in Ω) for s in SS)
    ES.S = Dict(s => [value(S[s,t]) for t in T] for s in SS)
    ES.F = Dict(l => [value(F[l,t]) for t in T] for l in Li)
    ES.F_adj = Dict(l => Dict(ω => [value(F_adj[l,t,ω]) for t in T] for ω in Ω) for l in Li)

    ES.B_up = Dict(p => Dict(ω => [value(B_up[p,t,ω]) for t in T] for ω in Ω) for p in P)
    ES.B_down = Dict(p => Dict(ω => [value(B_down[p,t,ω]) for t in T] for ω in Ω) for p in P);
    ES.B_up_L = Dict(s => Dict(ω => [value(B_up_L[s,t,ω]) for t in T] for ω in Ω) for s in SS)
    ES.B_down_L = Dict(s => Dict(ω => [value(B_down_L[s,t,ω]) for t in T] for ω in Ω) for s in SS);
    ES.B_up_G = Dict(s => Dict(ω => [value(B_up_G[s,t,ω]) for t in T] for ω in Ω) for s in SS)
    ES.B_down_G = Dict(s => Dict(ω => [value(B_down_G[s,t,ω]) for t in T] for ω in Ω) for s in SS);

    ES.λ_DA = Dict(n => [dual(NodeBalanceDA[n,t]) for t in T] for n in N)
    ES.λ_RT = Dict(n => Dict(ω => [dual(NodeBalanceRT[n,t,ω]) for t in T] for ω in Ω) for n in N)

    if ES.hydrogen_production_driver == "demand"
        ES.ρ_H2 = dual(HydrogenDemand)
    else
        ES.γ_H2 = sum(Electrolysers[e].η*(ES.L[e] .+ 
            1/length(Ω)*sum(ES.B_down[e][ω].-ES.B_up[e][ω] for ω in ES.Ω))
            for e in ES.E)
    end

    ES.model = m

end