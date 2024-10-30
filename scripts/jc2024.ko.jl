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
    met_ko = "HSST"    # Met
    ile_ko = "THRD_L"  # Ile
    
    # ko
    # bounds!(net0, met_ko, 0.0, 0.0)
    bounds!(net0, ile_ko, 0.0, 0.0)

    # med definition
    met_ex0 = "EX_met__L_e"
    ile_ex0 = "EX_ile__L_e"

    # Medium
    glc_ex0 = extras(net0, "EX_GLC")
    @show glc_ex0
    # bounds!(net0, glc_ex0, -10.0, 1000.0)
    # bounds!(net0, met_ex0, -1000.0, 1000.0)
    # bounds!(net0, met_ex0, -0.0, 1000.0)
    bounds!(net0, ile_ex0, -1000.0, 1000.0)
    # bounds!(net0, ile_ex0, -0.0, 1000.0)

    # objective functiom
    # max ( sum (ci * vi))
    # linear_weights!(net0, 0) # make all c = 0
    # I must be optimizing z1 + z2
    biom = extras(net0, "BIOM")
    linear_weights!(net0, biom, 1)
    
    # fba
    opm = FBAOpModel(net0, GLPK.Optimizer)
    optimize!(opm)
    @show objective_value(opm)
    @show solution(opm, biom)
    @show solution(opm, met_ex0)
    @show solution(opm, ile_ex0)
    @show solution(opm, glc_ex0)
end
