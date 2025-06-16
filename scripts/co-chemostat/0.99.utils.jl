## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
using MetX

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
# MARK: Double model
function double_model(net0, strain1, strain2, med = "Medium")
    
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
        ["$strain1:$rxn" for rxn in rxns0]; 
        ["$strain2:$rxn" for rxn in rxns0]; 
        ["$med:$rxn" for rxn in rxns0[ex_rxnis0]]; 
    ]
    @assert length(rxns1) == N1
    mets1 = [
        ["$strain1:$met" for met in mets0]; 
        ["$strain2:$met" for met in mets0]; 
        ["$med:$met" for met in mets0[ex_metis0]]; 
    ]
    @assert length(mets1) == M1
    b1 = zeros(M1)
    c1 = zeros(N1)

    # # Add biom auxiliar reaction
    # # TODO/ You can not just add specific growth rates
    # # z_pop = f_A z_A + f_B z_B
    # # - You can define a equal fraction growth
    # # - The populaation growth if both cells fractions are equal.

    # Θ2 = spzeros(M1, 1)
    # π0 = spzeros(1, N1)
    # for (rxni, rxn) in enumerate(rxns1)
    #     endswith(rxn, biom0) || continue
    #     π0[rxni] = 1
    # end
    # S2 = Float64[
    #     S1 Θ2 
    #     π0 -1
    # ]
    S2 = S1
    # lb2 = Float64[lb1; -1000.0]
    lb2 = lb1
    # ub2 = Float64[ub1; 1000.0]
    ub2 = ub1
    # rxns2 = [rxns1; "biom_aux_rxn"]
    rxns2 = rxns1
    # mets2 = [mets1; "biom_aux_met"]
    mets2 = mets1
    # b2 = [b1; 0]
    b2 = b1
    # c2 = [c1; 1]
    c2 = c1
    
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
    extras!(net2, "BIOM1", "$strain1:$biom0")
    extras!(net2, "BIOM2", "$strain2:$biom0")
    extras!(net2, "BIOM.AUX", "biom_aux_rxn")
    extras!(net2, "EX_RXNS0", rxns0[ex_rxnis0])
    extras!(net2, "EX_METS0", mets0[ex_metis0])
    
    extras!(net2, "is.double.net", true)
    extras!(net2, "strain1", strain1)
    extras!(net2, "strain2", strain2)
    extras!(net2, "med", med)

    return net2
end


# MARK: summary
function _opm_summary_PamModel_alter(opm) 
    # biom
    @show solution(opm, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f")
    @show solution(opm, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_b")
    @show solution(opm, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f")
    @show solution(opm, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_b")
    homo_pop_growth = 
        0.5 * solution(opm, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f") + 
        0.5 * solution(opm, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f")
    @show homo_pop_growth
    # glc
    println("-"^30)
    @show solution(opm, "ΔMet:EX_glc__D_e_f")
    @show solution(opm, "ΔMet:EX_glc__D_e_b")
    @show solution(opm, "ΔIle:EX_glc__D_e_f")
    @show solution(opm, "ΔIle:EX_glc__D_e_b")
    @show solution(opm, "Medium:EX_glc__D_e_f")
    @show solution(opm, "Medium:EX_glc__D_e_b")
    # met
    println("-"^30)
    @show solution(opm, "ΔMet:EX_met__L_e_f")
    @show solution(opm, "ΔMet:EX_met__L_e_b")
    @show solution(opm, "ΔIle:EX_met__L_e_f")
    @show solution(opm, "ΔIle:EX_met__L_e_b")
    @show solution(opm, "Medium:EX_met__L_e_f")
    @show solution(opm, "Medium:EX_met__L_e_b")
    # ile
    println("-"^30)
    @show solution(opm, "ΔMet:EX_ile__L_e_f")
    @show solution(opm, "ΔMet:EX_ile__L_e_b")
    @show solution(opm, "ΔIle:EX_ile__L_e_f")
    @show solution(opm, "ΔIle:EX_ile__L_e_b")
    @show solution(opm, "Medium:EX_ile__L_e_f")
    @show solution(opm, "Medium:EX_ile__L_e_b")

end

function _opm_summary_ecoli_core(opm) 
    # biom
    @show solution(opm, "BIOMASS_Ecoli_core_w_GAM")
    @show solution(opm, "EX_glc__D_e")
    @show solution(opm, "EX_nh4_e")
    @show solution(opm, "EX_lac__D_e")
    @show solution(opm, "EX_o2_e")
    @show solution(opm, "EX_co2_e")
end

# MARK: dict of vec (dov)
function dovGet(dict::Dict, id::String, dflt = nothing)
    return get(dict, id, dflt)
end

function dovPush!(dict::Dict, id::String, val) 
    vec = get!(dict, id, [])
    push!(vec, val)
    return val
end

function dovLast!(dict::Dict, id::String, dflt)
    vec = get!(dict, id, [])
    if isempty(vec)
        push!(vec, dflt)
    end
    return last(vec)
end

function dovLast(dict::Dict, id::String, dflt = nothing)
    vec = get!(dict, id, [])
    if isempty(vec)
        return dflt
    end
    return last(vec)
end

function _summaryNonZeroExchange(net::MetNet, opm::OpModel; th = 1e-2)
    # Get the non-zero values
    for rxn in reactions(net)
        contains(rxn, "EX_") || continue
        val = solution(opm, rxn)
        abs(val) < th && continue
        println("-"^30)
        summary(net, rxn)
        @show val
    end
end

# MARK: CLI
function _tryparse(str::AbstractString)
    for T in [Int, Float64, Bool]
        val = tryparse(T, str)
        isnothing(val) || return val
    end
    return str
end

CLI = Dict{String, Any}()
function _parse_ARGS(args=ARGS)
    for arg in args
        contains(arg, "=") || continue
        dig = split(arg, "=")
        length(dig) == 2 || continue
        CLI[dig[1]] = _tryparse(dig[2])
    end
    return
end

macro cli(ex::Expr)
    ex.head == Symbol("=") || 
        error("Expected assignation")
    sym = ex.args[1]
    symstr = string(sym)
    expr = ex.args[2]
    return quote
        $(esc(sym)) = get(CLI, $(esc(symstr)), $(esc(expr)))
    end
end


# This interface determine if an object is lite or not
# overwrite to change definition

# using Dates

# islite(::Any) = false    # fallback
# islite(::Number) = true
# islite(::Symbol) = true
# islite(::DateTime) = true
# islite(::VersionNumber) = true
# islite(::Nothing) = true
# islite(s::AbstractString) = length(s) < 1000

# # macro scope()
# #     return quote
# #         local _mod = @__MODULE__
# #         local _glob = Dict{Symbol, Any}()
# #         for f in names(_mod)
# #             # @show f
# #             # @show isdefined(_mod, f)
# #             _glob[f] = isdefined(_mod, f) ? getfield(_mod, f) : :UNDEFINED
# #         end
# #         local _loc = Base.@locals
# #         local _scope = Dict{Symbol, Any}()
# #         merge!(_scope, _glob, _loc)
# #         for (key, val) in _scope
# #             islite(val) && continue
# #             _scope[key] = hash(val)
# #         end
# #         _scope
# #     end
# # end


## -- .. - .-- .-. . .... -- -- -- .. ...
function _indexes(xs, ys;
        nsamples = 10_000,
        idx0 = 0.0,
        idx1 = 1.0,
    )
    idxs = intersect(eachindex(xs), eachindex(ys))
    nidxs = length(idxs)
    idx_min, idx_max = extrema(idxs)
    idx0 = clamp(floor(Int, idx0 * nidxs), idx_min, idx_max)
    idx1 = clamp(ceil(Int, idx1 * nidxs), idx_min, idx_max)
    idxs = idxs[idx0:idx1]
    if nsamples > 0
        return rand(idxs, nsamples)
    end
    return idxs
end

## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 

# TODO/ Move to Bloberias

import Base.rand
function Base.rand(rng::AbstractRNG, bb::BlobBatch)
    _uuids = getbuuids!(bb)
    _rid = rand(rng, _uuids)
    return blob(bb, _rid)
end

function Base.rand(bb::BlobBatch)
    rng = Random.default_rng()
    return rand(rng, bb)
end

function Base.rand(rng::AbstractRNG, B::Bloberia, bbid_prefix = nothing; )
    # eachbatch sort using all files from readdir
    bbs = eachbatch(B, bbid_prefix;
        sortfun = fns -> shuffle!(rng, fns),
        ch_size = 1,
        n_tasks = 1
    )
    return first(bbs)
end

function Base.rand(B::Bloberia, bbid_prefix = nothing)
    rng = Random.default_rng()
    return rand(rng, B, bbid_prefix)
end



## -- .- .- .-.-.- .-. - -.-. ..- .- -. -.. . 
nothing