# MARK: using
@time begin
    using EColi_Coexistence_2024
    using MetX
    using MetX.MetXNetHub
    using MetX.GLPK
    using MetX.Ipopt
    using MetX.Clp
    using Gurobi
    using CSV
    using DataFrames
    using JSON
    using CairoMakie
    using Serialization
    using SparseArrays
    using Random
end

# -- .. - .-- .-. . .... -- -- -- .. ...
include("0.utils.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
#=
DOING
- Reproduce @joyStudyGrowthEscherichia2010 results
- dfba simulations of a chemostat for a two strain system.
- testing different objective functions.
- testing different initial conditions.
=# 

# -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Prepare opm
include("1.1.prepare.PamModel-alter.2strains.opm.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: experimental data
let
    # json_file = joinpath(PaperSON_dir(), "raw.json")
    # global rawData = JSON.parsefile(json_file)
    # global fig1Data = rawData["data"]["Fig 1"]
    # nothing
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: DFBA
# dfba
SOL = Dict{String, Vector{Float64}}()
let
    cli = _parse_ARGS()

    # open medium, the glc will be controlled at each cell
    ub!(opm, "Medium:EX_glc__D_e_b", 1000.0)
    # open essential intake
    ub!(opm, "ΔMet:EX_ile__L_e_b", 1000.0)
    ub!(opm, "ΔIle:EX_ile__L_e_b", 1000.0)
    ub!(opm, "ΔMet:EX_met__L_e_b", 1000.0)
    ub!(opm, "ΔIle:EX_met__L_e_b", 1000.0)


    # MARK: ....objective function
    set_linear_obj!(opm, 
        [   
            "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_glc__D_e_b",  
            "ΔIle:EX_glc__D_e_b",
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ],
        [
              1e4,  # "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f" max ~1.0
              1e2,  # "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
              0.0,  # "Medium:EX_glc__D_e_b" 
            -1e-3,  # "ΔMet:EX_glc__D_e_b" max ~10.0
            -1e-3,  # "ΔIle:EX_glc__D_e_b"
            -1e-5,  # "ΔMet:EX_ile__L_e_b"
            -1e-5,  # "ΔIle:EX_ile__L_e_b"
            -1e-5,  # "ΔMet:EX_met__L_e_b"
            -1e-5,  # "ΔIle:EX_met__L_e_b"
        ]
    )
    
    # MARK: ....constants
    # from https://doi.org/10.3182/20100707-3-BE-2012.0059
    @cli glc_Km = 1.5               # [mM] 
    @cli glc_Vmax = abs(u_max)      # [mmol/g/h] Match PAM u_max
    @cli dt_min = 5e-3              # [h]
    @cli dt_max = 5e-2              # [h]
    @cli X_max = 10.0               # [g/L]
    @cli X_min = 1e-3               # [g/L]

    # MARK: ....initial cond
    @cli c_glc__D_e_t0 = 10.0
    dovPush!(SOL, "c.glc__D_e.ts", c_glc__D_e_t0)      # [mM]
    @cli s_glc__D_e_t0 = 10.0
    dovPush!(SOL, "s.glc__D_e.ts", s_glc__D_e_t0)       # [mM]
    # from https://doi.org/10.3182/20100707-3-BE-2012.0059 
    @cli ΔIle_X_t0 = 0.3
    dovPush!(SOL, "ΔIle:X.ts", ΔIle_X_t0)           # [g/L] 
    @cli ΔMet_X_t0 = 0.3
    dovPush!(SOL, "ΔMet:X.ts", ΔMet_X_t0)           # [g/L] 
    # TODO: make dt proportional to z
    @cli dt_t0 = 1e-1 # [h]
    dovPush!(SOL, "dt.ts", dt_t0)              # [h]
    @cli D_t0 = 0.2                         # [1/h]
    dovPush!(SOL, "D.ts", D_t0)                 # [1/h]
    
    @cli niters = 400
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
        # - 
        ub!(opm, "ΔIle:EX_glc__D_e_b", ΔIle_U_ex_glc_t) 
        ub!(opm, "ΔMet:EX_glc__D_e_b", ΔMet_U_ex_glc_t) 
        
        optimize!(opm)

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
        @show D_t
        @show sglc_t
        @show ΔIle_X_t
        @show ΔMet_X_t
        @show ΔIle_U_ex_glc_t
        @show ΔMet_U_ex_glc_t
        @show dt_t
        println("-"^30)
        _opm_summary_PamModel_alter(opm)

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

        # MARK: .... XA/XB adjust
        # I need to adjust all exchanges between the strains 
        # Given its X
        # EX_meti_b^cellB <= (EX_meti_f^cellA * X^cellA + c_meti D) / X^cellB
        # Here all c_meti = 0 mM
        
        # for bare_id in keys(carbon_sources)
        # bare_id == "EX_glc__D_e_f" && continue
        for bare_id in ["EX_ile__L_e", "EX_met__L_e"]
            @show bare_id
            ΔIle_ex_f_id = "ΔIle:$(bare_id)_f"
            ΔIle_ex_b_id = "ΔIle:$(bare_id)_b"
            ΔMet_ex_f_id = "ΔMet:$(bare_id)_f"
            ΔMet_ex_b_id = "ΔMet:$(bare_id)_b"

            ΔIle_ex_f = solution(opm, ΔIle_ex_f_id)
            ΔIle_ex_b = solution(opm, ΔIle_ex_b_id)
            ΔMet_ex_f = solution(opm, ΔMet_ex_f_id)
            ΔMet_ex_b = solution(opm, ΔMet_ex_b_id)

            _eps = 1e-2
            # if (ub(opm, ΔIle_ex_b_id) > 0)
                _U = (ΔMet_ex_f * ΔMet_X_t / ΔIle_X_t) + _eps
                @show _U
                ub!(opm, ΔIle_ex_b_id, _U)
            # end
            # if (ub(opm, ΔMet_ex_b_id) > 0)
                _U = (ΔIle_ex_f * ΔIle_X_t / ΔMet_X_t) + _eps
                @show _U
                ub!(opm, ΔMet_ex_b_id, _U)
            # end
        end
    end

    # Store
    datdir = string(@__FILE__, "data")
    mkpath(datdir)
    _scope = @scope
    filename = joinpath(datdir, string(hash(_scope), ".jls"))
    serialize(filename, Dict(
        "scope" => _scope,
        "SOL" => SOL,
        "opm" => opm, 
        "net0" => net0, 
    ))
end

