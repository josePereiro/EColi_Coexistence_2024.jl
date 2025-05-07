@time begin
    using MetX
    using MetX: MetXNetHub
    using GLPK
    using Ipopt
    using Gurobi
    using CairoMakie
    using SparseArrays
end

# TODO/TAI: remove obj function interface from the network

# -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
include("utils.jl")

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# ko
let
    global net0 = pull_net("iJO1366")
    bounds!(net0, "METabcpp", -1000.0, 1000.0)
    id1 = "Met"
    id2 = "Ile"
    med = "Medium"

    global net1 = double_model(net0, id1, id2; med)
    ko_id1 = "$id1:HSST"
    ko_id2 = "$id2:THRD_L"

    extras!(net1, "$(id1)_KO", ko_id1)
    extras!(net1, "$(id2)_KO", ko_id2)
    
    # close bounds
    bounds!(net1, ko_id1, 0.0, 0.0)
    bounds!(net1, ko_id2, 0.0, 0.0)

    # Open intra-cells exchanges
    exchs0 = extras(net1, "EX_RXNS0")
    for ex0 in exchs0
        bounds!(net1, "$id1:$ex0", -1000.0, 1000.0)
        bounds!(net1, "$id2:$ex0", -1000.0, 1000.0)
    end

    # some ids
    biom_aux = extras(net1, "BIOM.AUX")
    biom1 = extras(net1, "BIOM1")
    biom2 = extras(net1, "BIOM2")
    
    met_ex0 = "EX_met__L_e"
    met_id1 = "$id1:$met_ex0"
    met_id2 = "$id2:$met_ex0"
    met_med = "$med:$met_ex0"
    
    ile_ex0 = "EX_ile__L_e"
    ile_id1 = "$id1:$ile_ex0"
    ile_id2 = "$id2:$ile_ex0"
    ile_med = "$med:$ile_ex0"

    # Medium
    glc_ex0 = extras(net0, "EX_GLC")
    glc_id1 = "$id1:$glc_ex0"
    glc_id2 = "$id2:$glc_ex0"
    glc_med = "$med:$glc_ex0"
    @show glc_med
    bounds!(net1, glc_id1, -1000.0, 1000.0)
    bounds!(net1, glc_id2, -1000.0, 1000.0)
    bounds!(net1, glc_med, -10.0, 1000.0) # limiting

    bounds!(net1, met_id1, -1000.0, 1000.0)
    bounds!(net1, met_id2, -1000.0, 1000.0)
    bounds!(net1, met_med, - 0.0, 1000.0)  

    bounds!(net1, ile_id2, -1000.0, 1000.0)
    bounds!(net1, ile_id1, -1000.0, 1000.0)
    bounds!(net1, ile_med, -0.0, 1000.0) 

    # fba
    opm = FBAOpModel(lepmodel(net1), Ipopt.Optimizer)
    
    # max z1 + z2
    set_linear_obj!(opm, [biom1, biom2], [1, 1])
    optimize!(opm)
    tot_biom = objective_value(opm)
    
    # fix z1 + z2
    @show tot_biom
    lb!(opm, biom_aux, tot_biom - 1e-3)
    ub!(opm, biom_aux, tot_biom + 1e-3)

    global z1_vec = Float64[]
    global z2_vec = Float64[]
    global glc1_vec = Float64[]
    global glc2_vec = Float64[]

    global z1_fracv = range(0.0, 1.0; length = 5)
    for z1_frac in z1_fracv
        @show z1_frac
        # fix z1 fraction
        z1 = tot_biom * z1_frac
        bounds!(opm, biom1, z1 - 1e-3, z1 + 1e-3)

        # free glc
        bounds!(opm, glc_med, -1000.0, 1000.0) # limiting
        
        # min glc
        # set_linear_obj!(opm, [glc_id1, glc_id2], [1, 1])
        set_v2_obj!(opm, [glc_id1, glc_id2], MIN_SENSE)
        optimize!(opm)

        # record data
        push!(z1_vec, solution(opm, biom1))
        push!(z2_vec, solution(opm, biom2))
        push!(glc1_vec, solution(opm, glc_id1))
        push!(glc2_vec, solution(opm, glc_id2))
    end
    
    # min z1^2 + z2^2
    # set_v2_obj!(opm, [biom1, biom2], MIN_SENSE)
    
    # min glc^2 + glc^2
    # set_v2_obj!(opm, [biom1, biom2, glc_id1, glc_id2], MIN_SENSE)
    # optimize!(opm)

end

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
let
    # Solo esenciales
    glc1_vec = [-0.503243636915135, -2.776580115886399, -5.062104413993588, -7.347639979411194, -9.633152446986829]
    glc2_vec = [-9.486659294745063, -7.21380322976935, -4.928758784848706, -2.6437111366803205, -0.3604104323815521]

    # y intercept
    # -0.503243636915135 # Met
    # -0.3604104323815521 # Ile

    # biomass eq
    # ... + -0.290529 Ile + -0.153686 Met + ... => biomass
end
## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# Plots
let
    f = Figure()
    ax = Axis(f[1,1];
        title = "Intercambio abierto",
        # title = "Intercambio solo esenciales",
        limits = (0, 1, 0, 10), 
        xlabel = "z_M / (z_I + z_M)",
        ylabel = "glucosa (mmol/gDW h)",
    )

    lines!(ax, z1_fracv, -glc1_vec; color = :red, label = "ΔMetA")
    lines!(ax, z1_fracv, -glc2_vec; color = :blue,  label = "Δilv")
    axislegend(ax)
    f
end

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# #TODO: Make a package for Ploting with the sole proporse of avoiding
# update of ploting dependencies...
# Suffix: ex: Plots_LTS, CairoMakie_LTS

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 

    # # show
    # @show id1
    # @show id2
    # @show objective_value(opm)
    # @show solution(opm, biom1)
    # @show solution(opm, biom2)
    # @show solution(opm, met_id1)
    # @show solution(opm, met_id2)
    # @show solution(opm, met_med)
    # @show solution(opm, ile_id1)
    # @show solution(opm, ile_id2)
    # @show solution(opm, ile_med)
    # @show solution(opm, glc_med)
    # @show solution(opm, glc_id1)
    # @show solution(opm, glc_id2)

    # println("-."^20)
    # for exrxn0 in extras(net1, "EX_RXNS0")
    #     for id in [id1, id2]
    #         exrxn = "$id:$exrxn0"
    #         v = solution(opm, exrxn)
    #         abs(v) < 1e-3 && continue
    #         println(exrxn, "  ", v)
    #     end
    # end