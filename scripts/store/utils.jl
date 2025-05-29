## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# Double model
## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
function double_model(net0, id1, id2; med = "Medium")
    
    # S, lb, ub, rxns, mets = net0.S, net0.lb, net0.ub, net0.rxns, net0.mets
    # Network original fields
    S0 = stoi(net0)
    M0, N0 = size(S0)
    lb0, ub0 = bounds(net0)
    rxns0 = reactions(net0)
    mets0 = metabolites(net0)
    biom0 = extras(net0, "BIOM")

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

    # mix
    Θ0 = spzeros(M0, N0)
    E0 = S0[ex_metis0, ex_rxnis0]
    ME0, NE0 = size(E0)
    Θ1 = spzeros(M0, NE0)
    E1 = S0[ex_metis0, :] .* -1
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
    b1 = zeros(M1)
    c1 = zeros(N1)

    # Add biom auxiliar reaction
    Θ2 = spzeros(M1, 1)
    π0 = spzeros(1, N1)
    for (rxni, rxn) in enumerate(rxns1)
        endswith(rxn, biom0) || continue
        π0[rxni] = 1
    end
    S2 = Float64[
        S1 Θ2 
        π0 -1
    ]
    lb2 = Float64[lb1; -1000.0]
    ub2 = Float64[ub1; 1000.0]
    rxns2 = [rxns1; "biom_aux_rxn"]
    mets2 = [mets1; "biom_aux_met"]
    b2 = [b1; 0]
    c2 = [c1; 1]
    
    net2 = MetNet(;
        S=Matrix(S2), 
        lb=Vector(lb2), 
        ub=Vector(ub2), 
        rxns=rxns2, 
        mets=mets2,
        b=b2,
        c=c2
    )

    extras!(net2, "BIOM0", biom0)
    extras!(net2, "BIOM1", "$id1:$biom0")
    extras!(net2, "BIOM2", "$id2:$biom0")
    extras!(net2, "BIOM.AUX", "biom_aux_rxn")
    extras!(net2, "EX_RXNS0", rxns0[ex_rxnis0])
    extras!(net2, "EX_METS0", mets0[ex_metis0])
    
    extras!(net2, "is.double.net", true)
    extras!(net2, "id1", id1)
    extras!(net2, "id2", id2)
    extras!(net2, "med", med)

    return net2
end
