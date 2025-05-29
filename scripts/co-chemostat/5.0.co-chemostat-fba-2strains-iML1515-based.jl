## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: using
@time begin
    using MetX
    using MetX: MetXBase
    using MetX: MetXGEMs
    using MetX: MetXNetHub
    using MetX: Clp, GLPK, Ipopt
    using CSV
    using DataFrames
    using JSON
    using CairoMakie
    using Serialization
    using SparseArrays
    using Random
    using Scoperias
    using UUIDs
    using Dates
    using OrderedCollections
    using SparseArrays
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: include
include("0.0.project.jl")
include("0.88.scoperias.utils.jl")
include("0.99.utils.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("1.0.prepare.lepmodel.jl")
include("3.0.select.objfun.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
#=
DOING
Simulate a chemostat for two strains at different dilution rates
    See if it is the case that the dilution rate is not affecting coexistance
=# 


## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: sc_sel_hook!
sc_reset_sel_hooks!()
sc_sel_hook!([
    "co-chemostat-fba-2strains-iML1515-based", 
], 1) do scv::ScopeVariable

    startswith(scv.key, "_") && return :ignore
    return :include
end

sc_sel_hook!([
    "co-chemostat-fba-2strains-iML1515-based", 
    ["store.D"]
], 1) do scv::ScopeVariable
    startswith(scv.key, "_") && return :ignore
    return :include
end

# MARK: sc_call_hook!
_C = Dict{String, Any}() # RAM CACHE

sc_reset_call_hooks!()


# MARK: ..init_block
sc_call_hook!([
    "co-chemostat-fba-2strains-iML1515-based", 
    ["init_block"]
]) do sc::Scope

    empty!(_C)
    
    bb = sc["bb"][]
    rm(bb)

end

# MARK: ..prepare_lepmodel
sc_call_hook!([
    "co-chemostat-fba-2strains-iML1515-based", 
    "prepare_lepmodel"
]) do sc::Scope

    sc1 = filter(sc) do scv
        scv.key in [
            "blockver", "lepmodel_id", "lepmodel_args"
        ] && return true
        return false
    end
    
    # TODO/ make orderless hash function for scopes
    frame = repr(hash(sc1))
    @time _ref = blobio!(C, frame, "val", :getser!) do
        _net0, _lep0 = prepare_lepmodel(
            sc1["lepmodel_id"][], 
            sc1["lepmodel_args"][]
        )
        return (_net0, _lep0)
    end
    net0, lep0 = C[_ref]

    opm = get!(_C, frame) do
        FBAOpModel(lep0, Clp.Optimizer)
    end
    return (_ref, opm)
end

# MARK: ..store.D
sc_call_hook!([
    "co-chemostat-fba-2strains-iML1515-based", 
    ["store.D"]
]) do sc::Scope

    bb = sc["bb"][]
    b = blob!(bb, rbid())
    # TODO/ create an interface for labeling each variable to a destination
    # - similar to ScopeTape
    # - maybe reimplement ScopeTape but using Bloberias and Scoperias as deps

    merge!(b, "context", litecontext(sc))
    b["params", "lepmodel_args"] = sc["lepmodel_args"][]

    net0 = sc["net0"][]
    opm = sc["opm"][]

    sol = Dict{String, Float64}(
        rnx => solution(opm, rnx) 
        for rnx in reactions(net0)
    )
    b["sol", "."] = sol

end

# MARK: ..end
sc_call_hook!([
    "co-chemostat-fba-2strains-iML1515-based", 
    ["end"]
]) do sc::Scope

    bb = sc["bb"][]
    serialize!(bb)

    # finalizer
    empty!(bb)
    empty!(_C)
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Compute D serie
let
    # setup block
    @init_block "co-chemostat-fba-2strains-iML1515-based" "1"

    _rangeDs = range(0.1, 0.7; step = 1e-2)
    targetDs = rand(_rangeDs, 50)

    for (objfun_id, (lepmodel_id, lepmodel_args), D) in 
        Iterators.product(
            [
                "ΔIle.egoist.max.yield",
                "ΔMet.egoist.max.yield",
                "ΔIle.egoist.only",
                "ΔMet.egoist.only",
                "non-competitive.max.yield"
            ],
            [
                (
                    "iML1515.2strains", 
                    Dict()
                ),
                (
                    "PamModel-alter.2strains",
                    Dict(
                        "u_ider" => "EX_glc__D_e_b", 
                        "u_max" => 8.9
                    )
                )
            ],
            targetDs
        )

        prepare_lepmodel_ref, opm = @sc_call("prepare_lepmodel")
        net0, lep0 = C[prepare_lepmodel_ref]

        # SETUP
        objiders, objcoes = select_objfun(lepmodel_id, objfun_id)
        set_linear_obj!(opm, objiders, objcoes)

        _eps = 1e-3
        bounds!(opm, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f", D - _eps, D + _eps)
        bounds!(opm, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f", D - _eps, D + _eps)

        # OPTIMIZED
        try
            @time optimize!(opm)
        catch e
            (e isa InterruptException) && break
            continue
        end

        # INFO
        println("="^30)
        @show D
        @show objfun_id
        @show lepmodel_id
        @show objective_value(opm)
        
        _opm_summary_PamModel_alter(opm)

        println("-"^30)
        @show D / solution(opm, "ΔMet:EX_glc__D_e_b") 
        @show D / solution(opm, "ΔIle:EX_glc__D_e_b") 
        
        @sc_call "store.D"

    end

    @sc_call "end"

    nothing
end

