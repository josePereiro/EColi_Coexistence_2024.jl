@time begin
    using MetX
    using MetX: MetXNetHub
    using GLPK
    using SparseArrays
end

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# Double model
## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
function double_model(net0, id1, id2; med = "Medium")
    
    # S, lb, ub, rxns, mets = net0.S, net0.lb, net0.ub, net0.rxns, net0.mets
    # Network original fields
    global S0 = stoi(net0)
    global M0, N0 = size(S0)
    global lb0, ub0 = bounds(net0)
    global rxns0 = reactions(net0)
    global mets0 = metabolites(net0)

    # exchanges
    global ex_metis0 = Int[]
    global nonex_metis0 = Int[]
    for (meti, met) in enumerate(mets0)
        endswith(met, "_e") ? 
            push!(ex_metis0, meti) : 
            push!(nonex_metis0, meti)
    end

    global ex_rxnis0 = Int[]
    global nonex_rxnis0 = Int[]
    for (rxni, rxn) in enumerate(rxns0)
        startswith(rxn, "EX_") ? 
            push!(ex_rxnis0, rxni) : 
            push!(nonex_rxnis0, rxni) 
    end

    # mix
    Θ0 = spzeros(M0, N0)
    global E0 = S0[ex_metis0, ex_rxnis0]
    ME0, NE0 = size(E0)
    Θ1 = spzeros(M0, NE0)
    global E1 = S0[ex_metis0, :] .* -1
    E1[:,nonex_rxnis0] .= 0
    S1 = [
        S0 Θ0 Θ1
        Θ0 S0 Θ1
        E1 E1 E0
    ] |> dropzeros
    @assert sum(E0) == -sum(E1)
    M1, N1 = size(S1)
    # where:
    # S0: Original S0
    # Θ0: A matrix of zero of the proper size
    # Θ1: A matrix of zero of the proper size
    # E0: Submatrix of S0 with exchange reactions columns and the external metabolites rows
    # E1: Same information of E0 but exchange reactin are now located in the same indexes as in S0... Also, the coefficients are inverted

    # The rest is in correspondence with the definition of S1
    lb1 = [lb0; lb0; lb0[ex_rxnis0]]
    @assert length(lb1) == N1
    ub1 = [ub0; ub0; ub0[ex_rxnis0]]
    @assert length(ub1) == N1
    rxns1 = [
        ["$id1:$rxn" for rxn in rxns0]; 
        ["$id2:$rxn" for rxn in rxns0]; 
        ["$med:$rxn" for rxn in rxns0[ex_rxnis0]]; 
        
    ]
    @assert length(rxns1) == N1
    mets1 = [
        ["$id1:$met" for met in mets0]; 
        ["$id2:$met" for met in mets0]; 
        ["$med:$met" for met in mets0[ex_metis0]]; 
    ]
    @assert length(mets1) == M1
    
    global net1 = MetNet(;
        S=Matrix(S1), 
        lb=Vector(lb1), 
        ub=Vector(ub1), 
        rxns=rxns1, 
        mets=mets1,
        b=zeros(size(S1, 1)),
        c=zeros(size(S1, 2)),
    )

    biom0 = extras(net0, "BIOM")
    extras!(net1, "BIOM0", biom0)
    extras!(net1, "BIOM1", "$id1:$biom0")
    extras!(net1, "BIOM2", "$id2:$biom0")
    extras!(net1, "EX_RXNS0", rxns0[ex_rxnis0])
    extras!(net1, "EX_METS0", mets0[ex_metis0])
    
    extras!(net1, "is.double.net", true)
    extras!(net1, "id1", id1)
    extras!(net1, "id2", id2)
    extras!(net1, "med", med)

    return net1
end

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# ko
let
    global net0 = pull_net("iJO1366")
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

    # 
    exchs0 = extras(net1, "EX_RXNS0")
    for ex0 in exchs0
        bounds!(net1, "$id1:$ex0", -1000.0, 1000.0)
        bounds!(net1, "$id2:$ex0", -1000.0, 1000.0)
    end
    

    # med definition
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
    glc_id2 = "$id1:$glc_ex0"
    glc_med = "$med:$glc_ex0"
    @show glc_med
    bounds!(net1, glc_id1, -1000.0, 1000.0)
    bounds!(net1, glc_id2, -1000.0, 1000.0)
    bounds!(net1, glc_med, -10.0, 1000.0) # limiting

    bounds!(net1, met_id1, -1000.0, 1000.0)
    bounds!(net1, met_id2, -1000.0, 1000.0)
    bounds!(net1, met_med, -1000.0, 1000.0)  

    bounds!(net1, ile_id2, -1000.0, 1000.0)
    bounds!(net1, ile_id1, -1000.0, 1000.0)
    bounds!(net1, ile_med, -0.0, 1000.0) 

    # objective functiom
    # max ( sum (ci * vi))
    # linear_weights!(net1, 0) # make all c = 0
    # I must be optimizing z1 + z2
    biom1 = extras(net1, "BIOM1")
    biom2 = extras(net1, "BIOM2")
    linear_weights!(net1, 
        [biom1, biom2], 
        [0, 1]
    )
    
    # fba
    opm = FBAOpModel(lepmodel(net1), GLPK.Optimizer)
    optimize!(opm)
    @show id1
    @show id2
    @show objective_value(opm)
    @show solution(opm, biom1)
    @show solution(opm, biom2)
    @show solution(opm, met_id1)
    @show solution(opm, met_id2)
    @show solution(opm, met_med)
    @show solution(opm, ile_id1)
    @show solution(opm, ile_id2)
    @show solution(opm, ile_med)
    @show solution(opm, glc_med)
end
