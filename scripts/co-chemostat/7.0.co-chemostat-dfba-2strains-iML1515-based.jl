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
#=
DOING
- Reproduce @joyStudyGrowthEscherichia2010 results
- dfba simulations of a chemostat for a two strain system.
- testing different objective functions.
- testing different initial conditions.
=# 

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: include
include("0.0.project.jl")
include("0.88.scoperias.utils.jl")
include("0.99.utils.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("1.0.prepare.lepmodel.jl")
include("3.0.select.objfun.jl")

# -- .. - .-- .-. . .... -- -- -- .. ...
_parse_ARGS()

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: sc_sel_hook!
sc_reset_sel_hooks!()
sc_sel_hook!([
    "co-chemostat-dfba-2strains-iML1515-based", 
], 1) do scv::ScopeVariable

    startswith(scv.key, "_") && return :ignore
    return :include
end

sc_sel_hook!([
    "co-chemostat-dfba-2strains-iML1515-based", 
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
    "co-chemostat-dfba-2strains-iML1515-based", 
    ["init_block"]
]) do sc::Scope

    empty!(_C)
    
    bb = sc["bb"][]
    # rm(bb)

end

# MARK: ..prepare_lepmodel
sc_call_hook!([
    "co-chemostat-dfba-2strains-iML1515-based", 
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

# MARK: ..store.ts
sc_call_hook!([
    "co-chemostat-dfba-2strains-iML1515-based", 
    ["store.ts"]
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
    b["opm.last.sol", "."] = sol
    b["dyn.sol", "."] = deepcopy(sc["SOL"][])

    # DOING/ I need a rolling blobbatch
    # after a number of blobs I want a new one

    if length(bb) > 10
        blocklabel = sc["blocklabel"][]
        new_id = ttag_rname(blocklabel)
        @show new_id
        bb1 = rename_bb(bb, new_id)
        serialize!(bb1)
        empty!(bb)
    end

end

# MARK: ..end
sc_call_hook!([
    "co-chemostat-dfba-2strains-iML1515-based", 
    ["end"]
]) do sc::Scope

    bb = sc["bb"][]
    blocklabel = sc["blocklabel"][]
    new_id = ttag_rname(blocklabel)
    bb1 = rename_bb(bb, new_id)
    serialize!(bb1)

    # finalizer
    empty!(bb)
    empty!(_C)
end


## -- .. - .-- .-. . .... -- -- -- .. ...
# TODO/ create a simulation object so:
# - retain state
# - retain some record
# - allow restarting the simulation

# MARK: DFBA
# dfba
let
    # MARK: ....init_block
    @init_block "co-chemostat-dfba-2strains-iML1515-based" "1"
    
    while true

        # MARK: ....solution
        SOL = Dict{String, Vector{Float64}}()

        # MARK: ....setup model0            
        lepmodel_id, lepmodel_args = [
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
        ] |> rand
        
        prepare_lepmodel_ref, opm = @sc_call("prepare_lepmodel")
        net0, lep0 = C[prepare_lepmodel_ref]
        
        objfun_id = [
            "ΔIle.egoist.max.yield",
            "ΔMet.egoist.max.yield",
            # "ΔIle.egoist.only",
            # "ΔMet.egoist.only",
            "non-competitive.max.yield"
        ] |> rand
        objiders, objcoes = select_objfun(lepmodel_id, objfun_id)
        set_linear_obj!(opm, objiders, objcoes)
        
        # open medium, the glc will be controlled at each cell
        ub!(opm, "Medium:EX_glc__D_e_b", 1000.0)
        # open essential intake
        ub!(opm, "ΔMet:EX_ile__L_e_b", 1000.0)
        ub!(opm, "ΔIle:EX_ile__L_e_b", 1000.0)
        ub!(opm, "ΔMet:EX_met__L_e_b", 1000.0)
        ub!(opm, "ΔIle:EX_met__L_e_b", 1000.0)
        
        # MARK: ....constants
        # from https://doi.org/10.3182/20100707-3-BE-2012.0059
        u_max = [8.9] |> rand                # [mmol/g/h] Match PAM u_max
        glc_Km = [1.5] |> rand               # [mM] 
        glc_Vmax = [abs(u_max)] |> rand      # [mmol/g/h] Match PAM u_max
        dt_min, dt_max = [
            # (1e-4, 1e-2),
            (1e-2, 5e-2),
            # (5e-3, 5e-2),
        ] |> rand                            # [h]
        X_max = [10.0] |> rand               # [g/L]
        X_min = [1e-3] |> rand               # [g/L]

        # MARK: ....initial cond
        c_glc__D_e_t0 = [10.0] |> rand
        dovPush!(SOL, "c.glc__D_e.ts", c_glc__D_e_t0)      # [mM]
        s_glc__D_e_t0 = c_glc__D_e_t0
        dovPush!(SOL, "s.glc__D_e.ts", s_glc__D_e_t0)       # [mM]
        # from https://doi.org/10.3182/20100707-3-BE-2012.0059 
        X_t0_pool = Float64[ 0.3, 0.5, 0.8, 1.0 ]
        ΔIle_X_t0 = X_t0_pool |> rand
        dovPush!(SOL, "ΔIle:X.ts", ΔIle_X_t0)           # [g/L] 
        ΔMet_X_t0 = X_t0_pool |> rand
        dovPush!(SOL, "ΔMet:X.ts", ΔMet_X_t0)           # [g/L] 
        # TODO: make dt proportional to z
        dt_t0 = (dt_max + dt_min) * 0.5 # [h]
        dovPush!(SOL, "dt.ts", dt_t0)                   # [h]
        D_min = 0.1
        D_max = 0.6
        D_t0 = (D_max - D_min) * rand() + D_min         # [1/h]
        dovPush!(SOL, "D.ts", D_t0)                     # [1/h]
        abs_max_uglc_ref = nothing

        niters = 15_000
        for it in 1:niters

            # MARK: ....iter initial stuff
            ΔIle_X_t = dovLast(SOL, "ΔIle:X.ts")        # [g/L]
            ΔMet_X_t = dovLast(SOL, "ΔMet:X.ts")        # [g/L]
            sglc_t = dovLast(SOL, "s.glc__D_e.ts")      # [mM]
            cglc_t = dovLast(SOL, "c.glc__D_e.ts")      # [mM]
            D_t = dovLast(SOL, "D.ts")                  # [1/h]
            dt_t = dovLast(SOL, "dt.ts")                # [h]

            # MARK: ....break control
            if (sglc_t < 0) 
                println("sglc_t < 0!")
                break;
            end
            if (ΔIle_X_t < 0) 
                println("ΔIle_X_t < 0!")
                break;
            end
            if (ΔMet_X_t < 0) 
                println("ΔMet_X_t < 0!")
                break;
            end
            
            # MARK: ....compute bound(s)
            abs_max_uglc = (glc_Vmax * sglc_t) / (glc_Km + sglc_t)  # [mmol/g/h]
            # intake is a positive flux
            ΔIle_U_ex_glc_t = dovPush!(SOL, "ΔIle:U:EX_glc__D_e_b.ts", abs_max_uglc)   # [mmol/g/h] 
            ΔMet_U_ex_glc_t = dovPush!(SOL, "ΔMet:U:EX_glc__D_e_b.ts", abs_max_uglc)   # [mmol/g/h] 
            
            # MARK: ....fba
            # TAI/ Optimize only if U_ex_glc_t change significantly
            U_tol = 0.01
            do_optimize = 
                it < 5 || # head up
                abs(abs(abs_max_uglc_ref) - abs(abs_max_uglc)) / abs(abs_max_uglc) > U_tol 

            if (do_optimize)
                ub!(opm, "ΔIle:EX_glc__D_e_b", ΔIle_U_ex_glc_t) 
                ub!(opm, "ΔMet:EX_glc__D_e_b", ΔMet_U_ex_glc_t) 
                abs_max_uglc_ref = abs_max_uglc
                
                optimize!(opm)
            end

            # MARK: ....solutions
            # [mmol/g/h] 
            # *_v_ex_glc_t < 0 means intake
            ΔIle_v_ex_glc_t = 
                solution(opm, "ΔIle:EX_glc__D_e_f") - 
                solution(opm, "ΔIle:EX_glc__D_e_b")   
            dovPush!(SOL, "ΔIle:v:EX_glc__D_e_b.ts", ΔIle_v_ex_glc_t) 
            ΔMet_v_ex_glc_t = 
                solution(opm, "ΔMet:EX_glc__D_e_f") - 
                solution(opm, "ΔMet:EX_glc__D_e_b")   
            dovPush!(SOL, "ΔMet:v:EX_glc__D_e_b.ts", ΔMet_v_ex_glc_t) 
            
            # [1/h] 
            ΔIle_v_z_t = 
                solution(opm, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f") -
                solution(opm, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_b")
            dovPush!(SOL, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f", ΔIle_v_z_t) 
            ΔMet_v_z_t = 
                solution(opm, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f") -
                solution(opm, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_b")
            dovPush!(SOL, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f", ΔMet_v_z_t) 

            # MARK: ....info
            println(">"^30)
            if (do_optimize)
                println("optimize done!")
                @show it
                @show D_t
                @show sglc_t
                @show ΔIle_X_t
                @show ΔMet_X_t
                @show ΔIle_U_ex_glc_t
                @show ΔMet_U_ex_glc_t
                @show dt_t
                println("-"^30)
                _opm_summary_PamModel_alter(opm)
            else
                @show it
                @show abs_max_uglc
                @show sglc_t
            end

            # MARK: ....steps
            # X step
            ΔIle_X_vel = 0
            ΔIle_X_vel += ΔIle_v_z_t * ΔIle_X_t
            ΔIle_X_vel -= D_t * ΔIle_X_t
            ΔIle_X_tp1 = ΔIle_X_t + ΔIle_X_vel * dt_t # [g/L]
            ΔIle_X_tp1 = max(ΔIle_X_tp1, X_min)
            ΔMet_X_vel = 0
            ΔMet_X_vel += ΔMet_v_z_t * ΔMet_X_t
            ΔMet_X_vel -= D_t * ΔMet_X_t
            ΔMet_X_tp1 = ΔMet_X_t + ΔMet_X_vel * dt_t # [g/L]
            ΔMet_X_tp1 = max(ΔMet_X_tp1, X_min)
            
            # sglc
            sglc_vel = 0
            sglc_vel += cglc_t * D_t 
            sglc_vel -= sglc_t * D_t 
            sglc_vel += ΔIle_v_ex_glc_t * ΔIle_X_t 
            sglc_vel += ΔMet_v_ex_glc_t * ΔMet_X_t
            sglc_tp1 = sglc_t + sglc_vel * dt_t  # [mM]
            
            # update
            dovPush!(SOL, "ΔIle:X.ts", ΔIle_X_tp1)            
            dovPush!(SOL, "ΔMet:X.ts", ΔMet_X_tp1)            
            dovPush!(SOL, "s.glc__D_e.ts", sglc_tp1)            
            dovPush!(SOL, "c.glc__D_e.ts", cglc_t)   # constant
            dovPush!(SOL, "D.ts", D_t)                  # constant
            dt_tp1 = clamp(dt_max * sglc_tp1 / cglc_t, dt_min, dt_max)  # [h]
            dovPush!(SOL, "dt.ts", dt_tp1)

            # MARK: .... mass conservation
            # I need to adjust all exchanges between the strains 
            # Given its X to conserve mass.
            # EX_meti_b^cellB \le (EX_meti_f^cellA * X^cellA + c_meti D) / X^cellB
            # Here all c_meti = 0 mM
            
            # for bare_id in keys(carbon_sources)
            for bare_id in ["EX_ile__L_e", "EX_met__L_e"]
                # @show bare_id
                ΔIle_ex_f_id = "ΔIle:$(bare_id)_f"
                ΔIle_ex_b_id = "ΔIle:$(bare_id)_b"
                ΔMet_ex_f_id = "ΔMet:$(bare_id)_f"
                ΔMet_ex_b_id = "ΔMet:$(bare_id)_b"

                ΔIle_ex_f = solution(opm, ΔIle_ex_f_id)
                ΔIle_ex_b = solution(opm, ΔIle_ex_b_id)
                ΔMet_ex_f = solution(opm, ΔMet_ex_f_id)
                ΔMet_ex_b = solution(opm, ΔMet_ex_b_id)

                U_eps = 1e-2
                _U = (ΔMet_ex_f * ΔMet_X_t / ΔIle_X_t) + U_eps
                ub!(opm, ΔIle_ex_b_id, max(_U, 0.0))

                _U = (ΔIle_ex_f * ΔIle_X_t / ΔMet_X_t) + U_eps
                ub!(opm, ΔMet_ex_b_id, max(_U, 0.0))
            end
        
        end # for it in 1:niters
        
        @sc_call "store.ts"
        
    end # while true

    @sc_call "end"
end

