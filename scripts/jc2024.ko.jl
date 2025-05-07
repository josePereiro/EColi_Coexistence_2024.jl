@time begin
    using MetX
    using MetX: MetXNetHub
    using GLPK
    using SparseArrays
end

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# ko
let
    global net0 = pull_net("iJO1366")
    
    met_ex0 = "EX_met__L_e"
    ile_ex0 = "EX_ile__L_e"
    met_ko = "HSST"    # Met
    ile_ko = "THRD_L"  # Ile
    biom = extras(net0, "BIOM")
    
    # ko
    bounds!(net0, met_ko, 0.0, 0.0)
    # bounds!(net0, ile_ko, 0.0, 0.0)
    
    # Medium
    glc_ex0 = extras(net0, "EX_GLC")
    @show glc_ex0
    bounds!(net0, glc_ex0, -10.0, 1000.0)
    # bounds!(net0, met_ex0, -1000.0, 1000.0)
    bounds!(net0, met_ex0, -0.09, 1000.0)
    # bounds!(net0, ile_ex0, -0.1, 1000.0)
    # bounds!(net0, ile_ex0, -0.0, 1000.0)

    # objective functiom
    # max ( sum (ci * vi))
    # linear_weights!(net0, 0) # make all c = 0
    # I must be optimizing z1 + z2
    
    
    # fba
    opm = FBAOpModel(net0, GLPK.Optimizer)
    set_linear_obj!(opm, [biom], MAX_SENSE)
    optimize!(opm)
    @show objective_value(opm)
    @show solution(opm, biom)
    @show solution(opm, met_ex0)
    @show solution(opm, ile_ex0)
    @show solution(opm, glc_ex0)

    # dual price
    ex0 = met_ko
    test_points = lb(opm, ex0) .* [1.0, 0.95, 0.9]
    obj_m, obj_err, vars_ms, vars_errs = bound_dual_prices(
        opm, ex0, test_points, :lb; 
        dovars = true
    )
    obj_m

    # shadowprices
    # ile: -3.44199718444646
    # met: -6.497261404318072


    biom = [0.6497261404318072, 0.5847535263886264]
    ex0 = [-0.09999999999999998, -0.08999999999999997]
    (biom[2] - biom[1]) / (ex0[2] - ex0[1])
end
