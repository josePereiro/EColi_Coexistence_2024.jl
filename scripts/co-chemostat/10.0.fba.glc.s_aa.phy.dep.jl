## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
@time begin
    using CSV
    using MetX
    using MetX: Clp, GLPK, Ipopt
    using JSON
    using Random
    using MetXNetHub
    using DataFrames
    using CairoMakie
    using SparseArrays
    using Serialization
    using EColi_Coexistence_2024
    using Scoperias
    using Bloberias
    using Dates
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# DOING
# Computing the dependency of biomass with glc and aa concentrations

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: include
include("0.0.project.jl")
include("0.88.scoperias.utils.jl")
include("0.99.utils.jl")

# --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("1.0.prepare.lepmodel.jl")
include("3.0.select.objfun.jl")

# --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: sc_sel_hook!
sc_reset_sel_hooks!()
sc_sel_hook!([
    "glc.s_aa.phy", 
], 1) do scv::ScopeVariable

    startswith(scv.key, "_") && return :ignore
    return :include
end

# --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: sc_call_hook!

# MARK: ..store.pixel
sc_reset_call_hooks!()
sc_call_hook!([
    "glc.s_aa.phy", 
    ["store.pixel"]
]) do sc::Scope

    bb = sc["bb"][]
    b = blob!(bb, rbid())
    
    merge!(b, "context", litecontext(sc))
    b["params", "lepmodel_args"] = sc["lepmodel_args"][]

    net0 = sc["net0"][]
    opm = sc["opm"][]
    # TODO/ implement solution_dict
    sol = Dict{String, Float64}(
        rnx => solution(opm, rnx) 
        for rnx in reactions(net0)
    )
    b["opm.sol", "."] = sol

end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
michaelis(s, Vmax, KM) = Vmax * s / (KM + s)


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: Bound heatmap
let
    @init_block "fba.glc.q_aa.phy" "1"

    KM_glc = 1.75e-3 # mM
    KM_met = 0.13e-3 # mM
    KM_ile = 1.2e-3 # mM

    Vmax_glc = 8.9   # mmol/gh
    Vmax_met = 0.2   # mmol/gh
    Vmax_ile = 0.2   # mmol/gh

    rm(bb)
    
    SOL = Dict()

    glc_Umin, glc_Umax = 0.0, 10.0
    glc_Urange = range(glc_Umin, glc_Umax; length = 20)
    
    aa_Umin, aa_Umax = 0.0, 0.2
    aa_Urange = range(aa_Umin, aa_Umax; length = 20)

    for (aux_id, aux_ex, glc_ex) in [
        ("ΔIle", "EX_ile__L_e_b", "EX_glc__D_e_b")
        ("ΔMet", "EX_met__L_e_b", "EX_glc__D_e_b")
    ]
        lepmodel_args = Dict(
            "u_ider" => "EX_glc__D_e_b", 
            "u_max" => 8.9, 
            "aux" => aux_id,
        )
        net0, lep0 = prepare_lepmodel(
            "PamModel-alter.simple.aux", lepmodel_args
        )
        opm = FBAOpModel(lep0, Clp.Optimizer)

        for (glc_U, aa_U) in Iterators.product(glc_Urange, aa_Urange)
        
            bounds!(opm, aux_ex, 0.0, aa_U)
            bounds!(opm, glc_ex, 0.0, glc_U)

            # max z1 + z2
            set_linear_obj!(opm, 
                [
                    "BIOMASS_Ec_iML1515_WT_75p37M_f", 
                    "EX_glc__D_e_b",
                    "EX_met__L_e_b",
                    "EX_ile__L_e_b",
                ], 
                [1e4, -1e0, -1e-2, -1e-2]
            )
            @time optimize!(opm)

            # collect
            @sc_call "store.pixel"
        end

    end

    serialize!(bb)

    nothing
end