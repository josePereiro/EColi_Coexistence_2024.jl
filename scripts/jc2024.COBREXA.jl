## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
@time begin
    using GLPK
    using COBREXA
    const A = COBREXA.A
    using JSONFBCModels
    using LinearAlgebra
    using GLPK
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
# Load model
let
    fn = joinpath(@__DIR__, "iJO1366.json")
    global model = load_model(fn)
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
let
    # Network original fields
    global S0 = A.stoichiometry(model)
    global N0, M0 = size(S0)
    global lb0, ub0 = A.bounds(model)
    global rxns0 = A.reactions(model)
    global mets0 = A.metabolites(model)

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

    # exchangeless network
    global S1 = S0[nonex_metis0, nonex_rxnis0]
    global N1, M1 = size(S1)
    global lb1, ub1 = lb0[nonex_rxnis0], ub0[nonex_rxnis0]
    global rxns1 = rxns0[nonex_rxnis0]
    global mets1 = mets0[nonex_metis0]

    # mix
    Θ1 = zeros(N1, M1)
    E1 = S0[ex_metis0, ex_rxnis0]
    NE1, ME1 = size(E1)
    Θ2 = zeros(N1, ME1)
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
        ["A:$rxn" for rxn in rxns0[nonex_rxnis0]]; 
        ["B:$rxn" for rxn in rxns0[nonex_rxnis0]]; 
        rxns0[ex_rxnis0]
    ]
    @assert length(rxns2) == N2
    mets2 = [
        ["A:$met" for met in mets0[nonex_metis0]]; 
        ["B:$met" for met in mets0[nonex_metis0]]; 
        mets0[ex_metis0]
    ]
    @assert length(mets2) == M2
    

    # TODO: Brito: Como crear una red en COBREXA a partir de:
    # S2, lb2, ub2, rxns2, mets2

    S2

end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
let
    

    for met in mets
        endswith(met, "_e") || continue
        @show met
    end
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
# get stoichiometry matrix
let
    S = A.stoichiometry(model)
    A.reactions(model)
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
# get bounds
let
    lb, ub = A.bounds(model)
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 
# fba
let
    solution = flux_balance_analysis(model, optimizer = GLPK.Optimizer)
    solution[:fluxes]
end

## .- .-. - -.-- - - -. .. -.- .-.- .-. - --- - .- -.- - 