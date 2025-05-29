@time begin
    using MetX
    using MetX: MetXNetHub
    using GLPK
    using SparseArrays
end

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# Double model
## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
function double_model(net0, id1, id2)
    
    # S, lb, ub, rxns, mets = net0.S, net0.lb, net0.ub, net0.rxns, net0.mets
    # Network original fields
    S0 = stoi(net0)
    M0, N0 = size(S0)
    lb0, ub0 = bounds(net0)
    rxns0 = reactions(net0)
    mets0 = metabolites(net0)

    # exchanges
    ex_metis0 = Int[]
    nonex_metis0 = Int[]
    for (meti, met) in enumerate(mets0)
        endswith(met, "_e") ? 
            push!(ex_metis0, meti) : 
            push!(nonex_metis0, meti)
    end

    ex_rxnis0 = Int[]
    nonex_rxnis0 = Int[]
    for (rxni, rxn) in enumerate(rxns0)
        startswith(rxn, "EX_") ? 
            push!(ex_rxnis0, rxni) : 
            push!(nonex_rxnis0, rxni) 
    end

    # exchangeless network
    S1 = S0[nonex_metis0, nonex_rxnis0]
    N1, M1 = size(S1)
    lb1, ub1 = lb0[nonex_rxnis0], ub0[nonex_rxnis0]
    rxns1 = rxns0[nonex_rxnis0]
    mets1 = mets0[nonex_metis0]

    # mix
    Θ1 = spzeros(N1, M1)
    E1 = S0[ex_metis0, ex_rxnis0]
    NE1, ME1 = size(E1)
    Θ2 = spzeros(N1, ME1)
    E2 = S0[ex_metis0, nonex_rxnis0]
    S2 = [
        S1 Θ1 Θ2
        Θ1 S1 Θ2
        E2 E2 E1
    ]
    # where:
    # S1: Original S0 but removing exchange reactions columns and external metabolites rows
    # Θ1: A matrix of zero of the proper size
    # Θ2: A matrix of zero of the proper size
    # E1: Submatrix of S0 with exchange reactions columns and the external metabolites rows
    # E2: Submatrix of S0 with non exchange reactions columns and the external metabolites rows

    # The rest is in correspondence with the definition of S2
    M2, N2 = size(S2)
    lb2 = [lb0[nonex_rxnis0]; lb0[nonex_rxnis0]; lb0[ex_rxnis0]]
    @assert length(lb2) == N2
    ub2 = [ub0[nonex_rxnis0]; ub0[nonex_rxnis0]; ub0[ex_rxnis0]]
    @assert length(ub2) == N2
    rxns2 = [
        ["$id1:$rxn" for rxn in rxns0[nonex_rxnis0]]; 
        ["$id2:$rxn" for rxn in rxns0[nonex_rxnis0]]; 
        rxns0[ex_rxnis0]
    ]
    @assert length(rxns2) == N2
    mets2 = [
        ["$id1:$met" for met in mets0[nonex_metis0]]; 
        ["$id2:$met" for met in mets0[nonex_metis0]]; 
        mets0[ex_metis0]
    ]
    @assert length(mets2) == M2
    
    # TODO: Brito: Como crear una red en COBREXA a partir de:
    net2 = MetNet(;
        S=Matrix(S2), 
        lb=Vector(lb2), 
        ub=Vector(ub2), 
        rxns=rxns2, 
        mets=mets2,
        b=zeros(size(S2, 1)),
        c=zeros(size(S2, 2)),
    )

    biom0 = extras(net0, "BIOM")
    extras!(net2, "BIOM0", biom0)
    extras!(net2, "BIOM1", "$id1:$biom0")
    extras!(net2, "BIOM2", "$id2:$biom0")
    
    extras!(net2, "is.double.net", true)
    extras!(net2, "id1", id1)
    extras!(net2, "id2", id2)

    return net2
end

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# ko
let
    net0 = pull_net("iJO1366")
    id1 = "Met"
    id2 = "Ile"
    net2 = double_model(net0, id1, id2)
    ko_id1 = "$id1:HSST"
    ko_id2 = "$id2:THRD_L"

    met_ex = "EX_met__L_e"
    met_id1 = "$id1:$met_ex"
    met_id2 = "$id2:$met_ex"
    
    ile_ex = "EX_ile__L_e"
    met_id1 = "$id1:$ile_ex"
    met_id2 = "$id2:$ile_ex"

    extras!(net2, "$(id1)_KO", ko_id1)
    extras!(net2, "$(id2)_KO", ko_id2)
    
    # close bounds
    bounds!(net2, ko_id1, 0.0, 0.0)
    bounds!(net2, ko_id2, 0.0, 0.0)

    # medium definition
    glc_id0 = extras(net0, "EX_GLC")
    bounds!(net2, glc_id0, -10.0, 1000.0)

    # objective functiom
    # max ( sum (ci * vi))
    # linear_weights!(net2, 0) # make all c = 0
    # I must be optimizing z1 + z2
    biom1 = extras(net2, "BIOM1")
    biom2 = extras(net2, "BIOM2")
    linear_weights!(net2, 
        [biom1, biom2], 
        [1, 1]
    )
    
    # fba
    opm = FBAOpModel(lepmodel(net2), GLPK.Optimizer)
    optimize!(opm)
    @show id1
    @show id2
    @show objective_value(opm)
    @show solution(opm, biom1)
    @show solution(opm, biom2)
    @show solution(opm, met_id1)
    @show solution(opm, met_id2)
    @show solution(opm, ile_id1)
    @show solution(opm, ile_id2)
end
