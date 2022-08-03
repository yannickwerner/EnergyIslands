# function plot_electrolyser_operation(n, T, Ω, ES)

#     e = n*ifelse(occursin("HUB", n), "_ptg_offshore", "_ptg_onshore")
#     f = plot(
#         legend=:outerbottom,
#         title=e*" operation"
#         )

#     plot!(f, T, ES.L[e][T], label="Power demand DA")#, color=:black)
#     plot!(f, T, ES.L[e][T] .+ mean([ES.B_down[e][ω][T] -
#         ES.B_up[e][ω][T] for ω in Ω]), label="Expected power demand RT")
    
#     return f
# end


function plot_RES_generation(j, T, Ω, ES)
    fz=14
    f = plot(
        legend=false,
        xlabel="Timestep",
        ylabel="Power generation in GW",
        xtickfontsize=fz,
        ytickfontsize=fz,
        xguidefontsize=fz,
        yguidefontsize=fz,
        legendfontsize=fz,
    )
    for ω in Ω
        plot!(
            f,
            1:length(T),
            ES.Renewables[j].real[ω][T]./1000,
            label=ω,
            yguidefontsize=fz,)
    end
    return f
end


# function plot_capacity_factors(ES)
#     nodes = [ES.Electrolysers[e].node for e in ES.E]
#     capacity_factor_DA = [ES.Electrolysers[e].capacity_factors["DA"] for e in ES.E]
#     capacity_factor_BA_down = [ES.Electrolysers[e].capacity_factors["BA_down"] for e in ES.E]
#     capacity_factor_BA_up = [ES.Electrolysers[e].capacity_factors["BA_up"] for e in ES.E]

#     f = groupedbar([capacity_factor_BA_up capacity_factor_BA_down capacity_factor_DA],
#         bar_position=:stack,
#         bar_width=0.7,
#         xticks=(1:length(nodes), nodes),
#         xrotation = 45,
#         label=["BA_up" "BA_down" "DA"],
#         ylabel="Capacity factor",
#         title="Capacity factors electrolysers",
#         legend= :outerbottom
#         )

#     return f
# end


# function plot_capacity_factors_grouped(results)
#     nodes = []
#     val = AbstractFloat[]
#     market = []
#     group = []
#     market_dict = Dict(
#         "DA" => "Day-ahead dispatch",
#         "BA_up" => "Real-time upward adjustment",
#         "BA_down" => "Real-time downward adjustment",
#     )
#     for sc in keys(results)
#         ES = results[sc]
#         vars = ["DA", "BA_down", "BA_up"]
#         for var in vars
#             val = vcat(val, 
#                 [ES.Electrolysers[e].capacity_factors[var]
#                 for e in ES.E])
#             nodes = vcat(nodes, [ES.Electrolysers[e].node for e in ES.E])
#             market = vcat(market, repeat([market_dict[var]], length(ES.E)))
#             group = vcat(group, repeat([sc], length(ES.E)))
#         end
#     end

#     df = DataFrame(
#         :nodes => nodes,
#         :val => val,
#         :market => market,
#         :group => group
#     )

#     p = df |> @vlplot(
#         :bar,
#         x={:"group", axis={title=""}},
#         y={:"val", axis={title="Capacity factor"}},
#         color={
#             :"market",
#             axis={title=""},
#             legend={orient="bottom",
#                 direction="vertical"},
#         },
#         column={:"nodes", axis={title=""}},
#         spacing=2,
#         order="ordinal"
#         )
#         #width=800,height=800,
#     return p
# end


# function plot_capacity_factors_grouped_from_dataframe_pgf(df)
#     nodes = []
#     val = AbstractFloat[]
#     market = []
#     group = []
#     market_dict = Dict(
#         "DA" => "Day-ahead dispatch",
#         "BA_up" => "Real-time downward adjustment of electrolyser",
#         "BA_down" => "Real-time upward adjustment of electrolyser",
#     )

#     scenarios = unique(df[!, "scenario"])
#     for sc in scenarios
#         group_name = rsplit.(sc, "_")[3] #join(rsplit.(sc, "_")[1:3], "_")
#         # group_name = rsplit.(sc, "_")[6] #join(rsplit.(sc, "_")[1:3], "_")
#         countries = [n[1] for n in rsplit.(
#             df[df[!, :scenario] .== sc, :e], "_")]
#         no_e = length(countries)
#         vars = ["DA", "BA_down", "BA_up"]
#         for var in vars
#             val = vcat(val, 
#                 df[df[!, :scenario] .== sc, var])
#             nodes = vcat(nodes, countries)
#             market = vcat(market, repeat([market_dict[var]], no_e))
#             group = vcat(group, repeat([group_name], no_e))
#         end
#     end

#     df = DataFrame(
#         :nodes => nodes,
#         :val => val,
#         :market => market,
#         :group => group
#     )

#     p = df |> @vlplot(
#         :bar,
#         x={:"group", axis={title=""}},
#         y={:"val", axis={title="Capacity factor"}},
#         color={
#             :"market",
#             axis={title=""},
#             legend={orient="bottom",
#                 direction="horizontal",
#                 },
#         },
#         column={:"nodes", axis={title=""}},
#         spacing=5,
#         order="ordinal",
#         config={
#             axis={labelFontSize=20, titleFontSize=20},
#             legend={labelFontSize=20, titleFontSize=20},# "monospace", "titleFont": "monospace"},
#             # "header": {"labelFont": "monospace", "titleFont": "monospace"},
#             # "mark": {"font": "monospace"},
#             # "title": {"font": "monospace", "subtitleFont": "monospace"}
#             },
#         tickfontsize=22
#         )
#         #width=800,height=800,
#     return p
# end

# function plot_price_duration_curve(n, results)
#     f = plot(
#         xlabel="Hour",
#         ylabel="Day-ahead price in €/MWh"
#         )
#     for sc in keys(results)
#         price_sorted = sort(results[sc].λ_DA[n], rev=true)
#         plot!(f, price_sorted, label=sc)
#     end
#     return f
# end

# function plot_price_duration_curve_from_dataframe(n, result_dict)
#     f = plot(
#         # title=n,
#         xlabel="Hour",
#         ylabel="Day-ahead price in €/MWh"
#         )
#     for sc in keys(result_dict)
#         price_sorted = sort(result_dict[sc][!, "price_DA_$n"], rev=true)
#         label_sc = join(rsplit(sc, "_")[1:2], "_")
#         plot!(f, price_sorted, label=label_sc)
#     end
#     # plot!(f, sort(df[!, "price_DA_$n"], rev=true))
#     return f
# end



function plot_capacity_factors_grouped_from_dataframe(
    df,
    case)

    nodes = []
    val = AbstractFloat[]
    market = []
    group = []
    market_dict = Dict(
        "DA" => "Day-ahead dispatch electrolyser",
        "BA_up" => "Real-time downward adjustment electrolyser",
        "BA_down" => "Real-time upward adjustment electrolyser",
    )

    legend_orientation = "horizontal"

    bz_list = [
        "DKW1", "DKE1", "VindØ", "Bornholm", "NSWHP"]
    if case in ["BZ_comparison", "ntc_comparison"]
        flt = [n in bz_list for n in df[!,:n]]
        df = df[flt, :]
        legend_orientation = "vertical"
    end

    scenarios = unique(df[!, "scenario"])
    for sc in scenarios
        if case == "BZ_comparison"
            group_name = rsplit.(sc, "_")[1]
        elseif case == "year_comparison"
            group_name = rsplit.(sc, "_")[2]
        elseif case == "dataset_comparison"
            group_name = rsplit.(sc, "_")[3]
        elseif case == "ntc_comparison"
            group_name = rsplit.(sc, "_")[5]
        end

        countries = df[df[!, :scenario] .== sc, :n]
        no_e = length(countries)
        vars = ["DA", "BA_down", "BA_up"]
        for var in vars
            val = vcat(val, 
                df[df[!, :scenario] .== sc, var])
            nodes = vcat(nodes, countries)
            market = vcat(market, repeat([market_dict[var]], no_e))
            group = vcat(group, repeat([group_name], no_e))
        end
    end

    df = DataFrame(
        :nodes => nodes,
        :val => val,
        :market => market,
        :group => group
    )

    fz=12
    p = df |> @vlplot(
        :bar,
        x={
            :"group",
            axis={
                title="",
                labelFontSize=fz,
                titleFontSize=fz
            }
        },
        y={
            :"val",
            axis={
                title="Capacity factor",
                type="quantitative",
                labelFontSize=fz,
                titleFontSize=fz,
            }
        },
        color={
            :"market",
            axis={
                title="",              
                },
            legend={
                orient="bottom",
                direction=legend_orientation,
                labelFontSize=fz,
                titleFontSize=fz,
                labelLimit=500,
                },
        },
        column={
            :"nodes",
            axis={
                title="",
                labelFontSize=fz,
                titleFontSize=fz
            },
            sort="null"
        },
        spacing=2,
        order="ordinal"
        )
    return p
end

# function plot_price_duration_curve(n, results)
#     f = plot(
#         xlabel="Hour",
#         ylabel="Day-ahead price in €/MWh"
#         )
#     for sc in keys(results)
#         price_sorted = sort(results[sc].λ_DA[n], rev=true)
#         plot!(f, price_sorted, label=sc)
#     end
#     return f
# end

function plot_price_duration_curve_from_dataframe(n, result_dict)
    f = plot(
        xlabel="Hour",
        ylabel="Day-ahead price in €/MWh"
        )
    for sc in keys(result_dict)
        price_sorted = sort(result_dict[sc][!, "price_DA_$n"], rev=true)
        label_sc = join(rsplit(sc, "_")[1:2], "_")
        plot!(f, price_sorted, label=label_sc)
    end
    return f
end


# function plot_electrolyser_duration_curve(e, results)
#     f = plot(
#         xlabel="Hour",
#         ylabel="Power demand in MW"
#         )
#     for sc in keys(results)
#         price_sorted = sort(
#             results[sc].L[e] .+ mean([results[sc].B_down[e][ω] -
#             results[sc].B_up[e][ω] for ω in results[sc].Ω]), rev=true)
#         plot!(f, price_sorted, label=sc)
#     end
#     return f
# end


function plot_electrolyser_duration_curve_from_dataframe(e, result_dict)
    f = plot(
        xlabel="Hour",
        ylabel="Power demand in MW"
        )

    for sc in keys(result_dict)
        demand_sorted = sort(
            result_dict[sc][!, "L_$(e)"] .+
            result_dict[sc][!, "B_down_$(e)_avg"] .-
            result_dict[sc][!, "B_up_$(e)_avg"], rev=true)
        label_sc = join(rsplit(sc, "_")[1:2], "_")
        plot!(f, demand_sorted, label=label_sc)
    end
    return f
end


function plot_interconnector_flow_from_dataframe(I, ntc_dict, result_dict)

    fz=14
    f = plot(
        ylabel="Expected trade in MW",
        xlabel="Hour",
        xtickfontsize=fz,
        ytickfontsize=fz,
        xguidefontsize=fz,
        yguidefontsize=fz,
        legendfontsize=fz,
        grid=false,
        left_margin=2.5mm,
        xlim=[0,8760],
        )
    
    for sc in keys(result_dict)
        trade_dict = Dict(i => sort(
            result_dict[sc][!, "F_$(i)"] .+
            result_dict[sc][!, "F_adj_$(i)_avg"],
                rev=true)
            for i in I)
        sc_label = rsplit(sc, "_")[1]
        for i in I
            plot!(
                f,
                trade_dict[i],
                linewidth=2,
                label="$sc_label")
        end
    end

    years_unique = unique([rsplit(sc, "_")[2] for sc in keys(result_dict)])
    for y in years_unique
        for i in I
            # hline!(f, [ntc_dict[y][i]], label="Net-transfer capacity")
            plot!(
                f,
                [0, 8760],
                [ntc_dict[y][i], ntc_dict[y][i]],
                linewidth=2,
                label="Net-transfer capacity")
        end
    end

    return f
end

### FIX problem with length of horizon
function plot_price_during_congestion_from_dataframe(
    i, year, result_dict, ntc_dict; T=nothing)

    hub_rename_dict = Dict(
        "HUB1" => "VindØ",
        "HUB2" => "NSWHP",
        "HUB3" => "Bornholm"
    )

    fz=14
    f = plot(    
        legend=:bottomright,
        xlabel="Timestep",
        ylabel="Day-ahead price in €/MWh",
        xtickfontsize=fz,
        ytickfontsize=fz,
        xguidefontsize=fz,
        yguidefontsize=fz,
        legendfontsize=fz,
        grid=false,
        left_margin=2.5mm
        # xlim=[0,8760]
        )

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

            label = c in keys(hub_rename_dict) ?
                hub_rename_dict[c]*" "*rsplit(sc,"_")[1] :
                c*" "*rsplit(sc,"_")[1]

            plot!(
                f,
                result_dict[sc][!, "price_DA_$c"][hours_congestion][T],
                linewidth=2,
                label=label)
        end
    end
    return f
end


# function plot_nodal_system_balance(n, ω, T, ES)

#     value_spillage = sum([ES.G_spill[j][ω][T] for j in ES.J
#         if ES.Renewables[j].node == n])
#     label_spillage = "Spillage RES"    

#     export_lines = [l for l in ES.Li if ES.Interconnectors[l].from == n]
#     value_export = sum([ES.F_adj[l][ω][T] for l in export_lines])
#     import_lines = [l for l in ES.Li if ES.Interconnectors[l].to == n]
#     value_import = sum([-ES.F_adj[l][ω][T] for l in import_lines])

#     label_export = "Adjustment exports"
#     label_import = "Adjustment imports"

#     value_curtailment = ES.L_shed[n][ω][T]
#     label_curtailment = "Curtailment Load"
    
#     balance_values =
#         hcat(
#             value_spillage,
#             value_curtailment,
#             value_export,
#             value_import
#         )

#     balance_labels = 
#         hcat(
#             label_spillage,
#             label_curtailment,
#             label_export,
#             label_import
#         )

#     if ~isempty([e for e in ES.E if ES.Electrolysers[e].node == n])
#         value_electrolyser = sum([ES.B_down[e][ω][T] - ES.B_up[e][ω][T] for e in ES.E
#             if ES.Electrolysers[e].node == n])
#         label_electrolyser = "Balancing electrolyser"
#         balance_values = hcat(balance_values, value_electrolyser)
#         balance_labels = hcat(balance_labels, label_electrolyser)
#     end

#     if ~isempty([s for s in ES.SS if ES.Storages[s].node == n])
#         value_storage = sum([
#             ES.B_down_L[s][ω][T] + ES.B_up_L[s][ω][T] - (
#             ES.B_down_G[s][ω][T] + ES.B_up_G[s][ω][T])
#             for s in ES.SS 
#                 if ES.Storages[s].node == n])
#         label_storage = "Balancing storage" 
#         balance_values = hcat(balance_values, value_storage)
#         balance_labels = hcat(balance_labels, label_storage)
#     end

#     if ~isempty([r for r in ES.R if ES.Reservoirs[r].node == n])
#         value_reservoir = sum([
#             ES.B_down[r][ω][T] - ES.B_up[r][ω][T]
#             for r in ES.R
#                 if ES.Reservoirs[r].node == n]
#             )
#         label_reservoir = "Balancing reservoir"
#         balance_values = hcat(balance_values, value_reservoir)
#         balance_labels = hcat(balance_labels, label_reservoir)
#     end

#     if ~isempty([r for r in ES.R if ES.Reservoirs[r].node == n])
#         value_conventionals = sum([
#             ES.B_down[i][ω][T] - ES.B_up[i][ω][T]
#             for i in ES.I
#                 if ES.Conventionals[i].node == n])
#         label_conventionals = "Balancing conventionals"
#         balance_values = hcat(balance_values, value_conventionals)
#         balance_labels = hcat(balance_labels, label_conventionals)
#     end


#     f = groupedbar(
#         balance_values/1000,
#         bar_position=:stack,
#         bar_width=1,
#         label=balance_labels,
#         ylabel="Value in MW",
#         title="Balance in $n in scenario $ω",
#         legend=:outerbottom,
#         linecolor=:match,
#         grid=false,
#     )

#     deviation_res = 
#         sum([ES.Renewables[j].real[ω][T].-ES.G[j][T] for j in ES.J 
#             if ES.Renewables[j].node == n]) / 1000

#     plot!(
#         f,
#         deviation_res,
#         linewidth=1,
#         color=:black,
#         linetype=:stepmid,
#         label="RES deviation in scenario $ω"
#     )

#     return f
# end


function plot_nodal_system_balance_from_dataframe(n, ω, T, df)

    cols_res = rsplit.(
        names(df)[occursin.("G_spill_$n", names(df))],
            "_")
    names_res = unique([join(f[3:end-1], "_") for f in cols_res])
    value_spillage = sum([df[!, "G_spill_$(j)_$(ω)"] for j in names_res])
    label_spillage = "Spillage renewables"    

    export_lines = [name_split[2]*"_"*name_split[3] 
        for name_split in rsplit.(names(df[!, Not(:t)]), "_")
        if ((name_split[1] == "F") &
            (name_split[2] == n) &
            (length(name_split) == 3))]
    value_export = sum([df[!, "F_adj_$(i)_$(ω)"] for i in export_lines])
    label_export = "Adjustment exports"

    import_lines = [name_split[3]*"_"*name_split[2] 
    for name_split in rsplit.(names(df[!, Not(:t)]), "_")
    if ((name_split[1] == "F") &
        (name_split[2] == n) &
        (length(name_split) == 3))]
    value_import = sum([df[!, "F_adj_$(i)_$(ω)"] for i in import_lines])
    label_import = "Adjustment imports"

    value_curtailment = df[!, "L_shed_$(n)_$(ω)"]
    label_curtailment = "Curtailed load"
    
    balance_values =
        hcat(
            value_spillage,
            value_curtailment,
            value_export,
            value_import
        )

    balance_labels = 
        hcat(
            label_spillage,
            label_curtailment,
            label_export,
            label_import
        )

    if any(occursin.("L_$(n)_ptg", names(df)))
        e = join(rsplit.(
            names(df)[occursin.("L_$(n)_ptg",
            names(df))], "_")[1][2:4], "_")
        value_electrolyser =
            df[!, "B_down_$(e)_$(ω)"] .- df[!, "B_up_$(e)_$(ω)"]
        label_electrolyser = "Balancing electrolyser"
        balance_values = hcat(balance_values, value_electrolyser)
        balance_labels = hcat(balance_labels, label_electrolyser)
    end

    storages = []
    if any(occursin.("L_$(n)_pumped_storage", names(df))) 
        append!(storages, [join(rsplit.(
            names(df)[occursin.("L_$(n)_pumped_storage",
            names(df))], "_")[1][2:4], "_")])
    end
    if any(occursin.("L_$(n)_battery", names(df)))
        append!(storages, [join(rsplit.(
            names(df)[occursin.("L_$(n)_battery",
            names(df))], "_")[1][2:3], "_")]) 
    end
    if ~isempty(storages)
        value_storage = sum([
            df[!, "B_down_L_$(s)_$(ω)"] .+ df[!, "B_up_L_$(s)_$(ω)"] .-
            (df[!, "B_down_G_$(s)_$(ω)"] .+ df[!, "B_up_G_$(s)_$(ω)"])
            for s in storages])
        label_storage = "Balancing storage" 
        balance_values = hcat(balance_values, value_storage)
        balance_labels = hcat(balance_labels, label_storage)
    end

    if any(occursin.("G_$(n)_reservoir", names(df)))
        r = join(rsplit.(
            names(df)[occursin.("G_$(n)_reservoir",
            names(df))], "_")[1][2:3], "_")
        value_reservoir =
            df[!, "B_up_$(r)_$(ω)"] .- df[!, "B_down_$(r)_$(ω)"]
            label_reservoir = "Balancing reservoir"
            balance_values = hcat(balance_values, value_reservoir)
            balance_labels = hcat(balance_labels, label_reservoir)
    end



    names_conv = [conv for conv in names(df) if 
        occursin(n, conv) &
        (occursin("B_down_$(n)", conv) | occursin("B_up_$(n)", conv)) &
        occursin(ω, conv) &
        ~occursin("ptg", conv) &
        ~occursin("reservoir", conv)
        ]

    names_conv_down = names_conv[occursin.("down", names_conv)]
    names_conv_up = names_conv[occursin.("up", names_conv)]

    if ~isempty(names_conv)
        value_conventionals = sum([
            df[!, i] .- df[!, j]
            for (i,j) in zip(names_conv_up, names_conv_down)])
        label_conventionals = "Balancing conventionals"
        balance_values = hcat(balance_values, value_conventionals)
        balance_labels = hcat(balance_labels, label_conventionals)
    end

    fz=14
    f = groupedbar(
        balance_values[T,:],
        bar_position=:stack,
        bar_width=1,
        label=balance_labels,
        ylabel="Value in MW",
        title="Balance in $n in scenario $ω",
        legend=:outerbottom,
        linecolor=:match,
        grid=false,
        xtickfontsize=fz,
        ytickfontsize=fz,
        xguidefontsize=fz,
        yguidefontsize=fz,
        legendfontsize=fz,
    )

    deviation_res = sum([df[!, "Deviation_$(j)_$(ω)"][T] 
        for j in names_res])

    plot!(
        f,
        deviation_res,
        linewidth=1,
        color=:black,
        linetype=:stepmid,
        label="Forecast error renewables underlying scenario"
    )

    return f
end


# function plot_system_balance(ω, T, ES)

#     value_spillage = sum([ES.G_spill[j][ω][T] for j in ES.J])
#     label_spillage = "Spillage RES"    

#     value_curtailment = sum([ES.L_shed[n][ω][T] for n in ES.N])
#     label_curtailment = "Curtailment load"    

#     value_electrolyser = sum([ES.B_down[e][ω][T] - ES.B_up[e][ω][T] for e in ES.E])
#     label_electrolyser = "Balancing electrolyser"

#     value_storage = sum([
#         ES.B_down_L[s][ω][T] + ES.B_up_L[s][ω][T] - (
#         ES.B_down_G[s][ω][T] + ES.B_up_G[s][ω][T])
#         for s in ES.SS])
#     label_storage = "Balancing storage"

#     value_reservoir = sum([
#         ES.B_down[r][ω][T] - ES.B_up[r][ω][T]
#         for r in ES.R])
#     label_reservoir = "Balancing reservoir"

#     value_conventionals = sum([
#         ES.B_down[i][ω][T] - ES.B_up[i][ω][T]
#         for i in ES.I])
#     label_conventionals = "Balancing conventionals"

#     balance_values =
#         hcat(
#             value_spillage,
#             value_electrolyser,
#             value_curtailment,
#             value_storage,
#             value_reservoir,
#             value_conventionals
#         )./1000
    
#     balance_labels = 
#         hcat(
#             label_spillage,
#             label_electrolyser,
#             label_curtailment,
#             label_storage,
#             label_reservoir,
#             label_conventionals
#         )

#     fz = 14
#     fontsize_dict = Dict(
#         :xtickfontsize => fz,
#         :ytickfontsize => fz,
#         :xguidefontsize => fz,
#         :yguidefontsize => fz,
#         :legendfontsize => fz
#     )

#     f = groupedbar(
#         balance_values,
#         bar_position=:stack,
#         bar_width=1,
#         label=balance_labels,
#         ylabel="Value in GW",
#         xlabel="Hour",
#         title="System-wide balance in scenario $ω",
#         legend=:outerright,
#         linecolor=:match,
#         size=(1000,500),
#         bottom_margin = 7.5mm,
#         grid=false,
#         fontsize_dict...
#     )

#     deviation_res = 
#         sum([ES.Renewables[j].real[ω][T].-ES.G[j][T] for j in ES.J]) ./ 1000

#     plot!(
#         f,
#         deviation_res,
#         linewidth=1,
#         color=:black,
#         linetype=:stepmid,
#         label="RES deviation in scenario $ω"
#     )

#     return f
# end



function plot_system_balance_from_dataframe(ω, T, df)

    value_spillage = df[!, "system_spillage_RES_$ω"]
    label_spillage = "Spillage renewables"    

    value_curtailment = df[!, "system_curtailment_load_$ω"]
    label_curtailment = "Curtailed load"    

    value_electrolyser = df[!, "system_balancing_electrolyser_$ω"]
    label_electrolyser = "Balancing electrolyser"

    value_storage = df[!, "system_balancing_storage_$ω"]
    label_storage = "Balancing storage"

    value_reservoir = df[!, "system_balancing_reservoir_$ω"]
    label_reservoir = "Balancing reservoir"

    value_conventionals = df[!, "system_balancing_conventionals_$ω"]
    label_conventionals = "Balancing conventionals"

    balance_values =
        hcat(
            value_spillage,
            value_electrolyser,
            value_curtailment,
            value_storage,
            value_reservoir,
            value_conventionals
        )./1000
    
    balance_labels = 
        hcat(
            label_spillage,
            label_electrolyser,
            label_curtailment,
            label_storage,
            label_reservoir,
            label_conventionals
        )

    fz = 14

    f = groupedbar(
        balance_values[T,:],
        bar_position=:stack,
        bar_width=1,
        label=balance_labels,
        ylabel="Value in GW",
        xlabel="Timestep",
        xlim=[0,72],
        ylim=[-52,45],
        legend=:bottomleft,
        linecolor=:match,
        size=(1000,500),
        bottom_margin=5mm,
        left_margin=5mm,
        grid=false,
        xtickfontsize=fz,
        ytickfontsize=fz,
        xguidefontsize=fz,
        yguidefontsize=fz,
        legendfontsize=12
    )

    deviation_res = 
        df[!, "system_deviation_RES_$ω"][T] ./ 1000

    plot!(
        f,
        deviation_res,
        linewidth=1,
        color=:black,
        linetype=:stepmid,
        label="Forecast error renewables"
    )

    return f
end


# function plot_interconnector_congestion(I, T, ES)
#     trade_dict = Dict(i => sort(
#         ES.Interconnectors[i].expected_trade[T], rev=true)
#         for i in I)

#     fz = 12
#     f = plot(
#         ylabel="Expected trade in MW",
#         xlabel="Hour",
#         xtickfontsize=fz,
#         ytickfontsize=fz,
#         xguidefontsize=fz,
#         yguidefontsize=fz,
#         legendfontsize=fz
#         )
#     for i in I
#         ntc = ifelse(
#             i in ["DKW1_HUB1", "HUB1_DKW1"],
#             2000,
#             ES.Interconnectors[i].ntc
#         )
#         plot!(f, trade_dict[i], label=i)
#         hline!(f, [ntc], label="Capacity $i" )
#     end
#     return f
# end



# function plot_price_during_congestion(i, results; T=nothing)
#     hours_congestion = results["HBZ"].Interconnectors[i].expected_trade .>= 2000
#     if isnothing(T)
#         T = hours_congestion
#     end
#     f = plot(    
#         legend=:bottomright,
#         xlabel="Congested hours",
#         ylabel="Day-ahead price in €/MWh"
#         )        
#     for bz_config in ["HBZ", "OBZ"]
#         for c in ["DKW1", "HUB1"]
#             plot!(
#                 f,
#                 results[bz_config].λ_DA[c][hours_congestion][T],
#                 label=c*"_"*bz_config)
#         end
#     end

#     return f
# end

# f = plot(ES.T, ES.L["DELU_ptg_onshore"], label="DELU ptg", left_margin = left_margin, right_margin = right_margin)
# plot!(ES.T, ES.L["HUB1_ptg_offshore"], label="Hub1 ptg", color = :purple, left_margin = left_margin, right_margin = right_margin)
# plot!(twinx(), ES.T, ES.λ_DA["DELU"], legend=:bottomright, color = :red, label="Preis", left_margin = left_margin, right_margin = right_margin)
# # for ω in ES.Ω
# #     plot!(f, ES.T, ES.B_up["DELU_ptg_onshore"][ω], label=ω)
# # end
# @show f


# function plot_expected_nodal_balance(n, T, ES)
#     export_lines = [l 
#         for l in ES.Li if ES.Interconnectors[l].from == n]
#     import_lines = [l 
#         for l in ES.Li if ES.Interconnectors[l].to == n]

#     e = n*ifelse(occursin("HUB", n), "_ptg_offshore", "_ptg_onshore")
#     value_electrolyser = ES.L[e][T] .+ mean([ES.B_down[e][ω][T] -
#         ES.B_up[e][ω][T] for ω in ES.Ω])
#     label_electrolyser = "Expected electrolyser consumption"

#     j = n*"_offshore_wind"
#     value_windfarm = mean(ES.Renewables[j].real[ω][T] for ω in ES.Ω)
#     label_windfarm = "Expected windfarm production"

#     value_spillage = -mean(ES.G_spill[j][ω][T] for ω in ES.Ω)
#     label_spillage = "Expected windfarm spillage"

#     values = reduce(hcat,
#         vcat(
#             [value_windfarm],
#             [.-value_electrolyser],
#             [ES.Interconnectors[l].expected_trade[T] 
#                 for l in import_lines],
#             [-ES.Interconnectors[l].expected_trade[T]
#                 for l in export_lines],
#             [value_spillage]
#             ))

#     labels = reduce(hcat, vcat(
#         [label_windfarm],
#         [label_electrolyser],
#         "Expected trade " .* import_lines,
#         "Expected trade " .* export_lines,
#         [label_spillage]
#         ))

#     f = groupedbar(values,
#         bar_position=:stack,
#         bar_width=1,
#         label=labels,
#         ylabel="Value in MW",
#         title="Expected nodal balance $n",
#         legend= :outerbottom,
#         linecolor=:match,
#         grid=false,
#         )

#     return f
# end


# function plot_real_time_adjustment(n, T, ω, ES)

#     j = n*"_offshore_wind"

#     value_spillage = ES.G_spill[j][ω][T]
#     label_windfarm = "Spillage windfarm"

#     export_lines = [l for l in ES.Li if ES.Interconnectors[l].from == n]
#     value_export = [ES.F_adj[l][ω][T] for l in export_lines]
#     import_lines = [l for l in ES.Li if ES.Interconnectors[l].to == n]
#     value_import = [-ES.F_adj[l][ω][T] for l in import_lines]

#     balance_values = vcat([value_spillage],
#         value_import,
#         value_export
#         )
#     balance_labels = vcat(label_windfarm,
#         "Adjustment trade " .* import_lines,
#         "Adjustment trade " .* export_lines)

#     if ~(n in ["HUB1", "HUB2", "HUB3", "DKKF", "DEKF"])
#         value_curtailment = ES.L_shed[n][ω][T]
#         label_demand = "Curtailment demand"
#         balance_values = vcat(balance_values, [value_curtailment])
#         balance_labels = vcat(balance_labels, [label_demand])
#     end

#     if ~isempty([ES.Electrolysers[e].node == n for e in ES.E])
#         e = n*ifelse(occursin("HUB", n), "_ptg_offshore", "_ptg_onshore")
#         value_electrolyser = ES.B_down[e][ω][T] - ES.B_up[e][ω][T]
#         label_electrolyser = "Balancing electrolyser"
#         balance_values = vcat(balance_values, [value_electrolyser])
#         balance_labels = vcat(balance_labels, [label_electrolyser])
#     end

#     if ~isempty([ES.Storages[s].node == n for s in ES.SS])
#         value_storage = sum([
#             ES.B_down_L[s][ω][T] + ES.B_up_L[s][ω][T] - (
#             ES.B_down_G[s][ω][T] + ES.B_up_G[s][ω][T])
#             for s in ES.SS if ES.Storages[s].node == n])
#         label_storage = "Balancing storage"
#         balance_values = vcat(balance_values, [value_storage])
#         balance_labels = vcat(balance_labels, [label_storage])
#     end

#     if ~isempty([ES.Reservoirs[r].node == n for r in ES.R])
#         value_reservoir = [
#             ES.B_down[r][ω][T] - ES.B_up[r][ω][T]
#             for r in ES.R if ES.Reservoirs[r].node == n]
#         label_reservoir = "Balancing reservoir"
#         balance_values = vcat(balance_values, [value_reservoir])
#         balance_labels = vcat(balance_labels, [label_reservoir])
#     end

#     if ~isempty([ES.Conventionals[i].node == n for i in ES.I])
#         value_conventionals = [
#             ES.B_down[i][ω][T] - ES.B_up[i][ω][T]
#             for i in ES.I if ES.Conventionals[i].node == n]
#         label_conventionals = "Balancing conventionals"
#         balance_values = vcat(balance_values, [value_conventionals])
#         balance_labels = vcat(balance_labels, [label_conventionals])
#     end

#     f = groupedbar(
#         reduce(hcat, vcat(balance_values)),
#         bar_position=:stack,
#         bar_width=1,
#         label=reduce(hcat, vcat(balance_labels)),
#         ylabel="Value in MW",
#         title="Scenario specific balancing "*n,
#         legend=:outerbottom,
#         linecolor=:match,
#         grid=false,
#     )

#     plot!(
#         f,
#         ES.Renewables[j].real[ω][T].-ES.G[j][T],
#         linewidth=1,
#         color=:black,
#         linetype=:stepmid,
#         label="Wind deviation in scenario $ω"
#     )

#     return f
# end

# function plot_merit_order(ES)
#     mat = permutedims(hcat([[c, ES.Conventionals[c].g_max, ES.Conventionals[c].mc]
#         for c in ES.I]...))
#     df = DataFrame(:cap => mat[:,2], :mc => mat[:,3], :name => mat[:,1])
#     sort!(df, [:mc])
#     df[!, :cap_cum] = cumsum(df[!,:cap])./1000
#     f = plot(df[!, :cap_cum], df[!, :mc])
#     return f
# end


# function create_plots(n, ω, ES, T, Ω, path="./results/plots")
    
#     time = Dates.format(Dates.now(), "yyyy-mm-dd-HH-MM")
#     # config = "$(ES.year_model)_$(ES.bidding_zone_config)_$time"
#     config = "$(ES.year_model)_HBZ_$time"

#     f = plot_electrolyser_operation(n, T, Ω, ES)
#     savefig(f, path*"electrolyser_operation_$config.png")
#     f = plot_capacity_factors(ES)
#     savefig(f, path*"capacity_factors_$config.png")
#     f = plot_expected_nodal_balance(n, T, ES)
#     savefig(f, path*"expected_nodal_balance_$config.png")
#     f = plot_real_time_adjustment(n, T, ω, ES)
#     savefig(f, path*"real_time_adjustment_$config.png")
# end

