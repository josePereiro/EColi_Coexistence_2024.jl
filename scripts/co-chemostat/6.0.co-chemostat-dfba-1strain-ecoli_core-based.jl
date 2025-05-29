@time begin
    using EColi_Coexistence_2024
    using MetX
    using MetXNetHub
    using GLPK
    using Ipopt
    using CSV
    using DataFrames
    using JSON
    using Serialization
    using SparseArrays
    using Random
end

# -- .. - .-- .-. . .... -- -- -- .. ...
#=
DOING
- Reproduce @joyStudyGrowthEscherichia2010 results
- dfba simulations of a chemostat for a two strain system.
- testing different objective functions.
- testing different initial conditions.
=# 

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("0.0.project.jl")
include("0.99.utils.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("1.0.prepare.lepmodel.jl")
include("3.0.select.objfun.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: tape hooks
st_config!(r"co-chemostat-dfba-1strain-ecoli_core",
    "do.write.batch.at.len" => 300
)

st_hook!(r"co-chemostat-dfba-1strain-ecoli_core") do scv, sc

    ts_isclass(scv, :unknown) || return :pass

    # isa(scv.val, AbstractArray) && return :blob
    startswith(scv.key, "net") && return :blob
    startswith(scv.key, "opm") && return :blob

    scv.key == "args" && return :blob

    return :ignored
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: DFBA
# dfba
SOL = Dict{String, Vector{Float64}}()
let
    @st_flush! "co-chemostat-dfba-1strain-ecoli_core"

    lepmodel_id = "core_ecoli"
    net0, opm = prepare_lepmodel(lepmodel_id)


    # objective function
    objfun_id = "max.biomass.max.yield"
    objiders, objcoes = select_objfun(lepmodel_id, objfun_id)
    set_linear_obj!(opm, objiders, objcoes)
    
    # constants
    # from https://doi.org/10.3182/20100707-3-BE-2012.0059
    glc_Km = 0.015      # [mM] 
    glc_Vmax = 10.0     # [mmol/g/h]
    dt_max = 1e-1       # [h]
    dt_min = 1e-4       # [h]
    
    # initial cond
    dovPush!(SOL, "s.EX_glc__D_e.ts", 10.0)     # [mM]
    dovPush!(SOL, "c.EX_glc__D_e.ts", 10.0)     # [mM]
    # from https://doi.org/10.3182/20100707-3-BE-2012.0059
    dovPush!(SOL, "X.ts", 0.003)                # [g/L] 
    # TODO: make dt proportional to z
    dovPush!(SOL, "dt.ts", dt_min)              # [h]
    dovPush!(SOL, "D.ts", 3e-1)                 # [1/h]
    
    # "EX_nh4_e"
    # "EX_o2_e"
    # "EX_co2_e"
    
    niters = 100000
    for it in 1:niters
        # iter initial stuff
        X_t = dovLast(SOL, "X.ts")                   # [g/L]
        sglc_t = dovLast(SOL, "s.EX_glc__D_e.ts")    # [mM]
        cglc_t = dovLast(SOL, "c.EX_glc__D_e.ts")    # [mM]
        D_t = dovLast(SOL, "D.ts")                   # [1/h]
        dt_t = dovLast(SOL, "dt.ts")                 # [h]

        # break
        if (sglc_t < 0) 
            println("sglc_t < 0!")
            break;
        end
        
        # compute bound
        abs_max_uglc = (glc_Vmax * sglc_t) / (glc_Km + sglc_t)  # [mmol/g/h]
        # intake is a negative flux
        max_uglc_t = dovPush!(SOL, "max_u.EX_glc__D_e.ts", -abs_max_uglc)   # [mmol/g/h] 
        
        # fba
        lb!(opm, "EX_glc__D_e", max_uglc_t) 
        optimize!(opm)

        v_ex_glc_t = solution(opm, "EX_glc__D_e")       # [mmol/g/h] 
        dovPush!(SOL, "v.EX_glc__D_e.ts", v_ex_glc_t) 
        
        v_z_t = solution(opm, "BIOMASS_Ecoli_core_w_GAM")   # [1/h] 
        dovPush!(SOL, "v.BIOMASS_Ecoli_core_w_GAM.ts", v_z_t) 

        # println("-"^30)
        # @show sglc_t
        # @show v_ex_glc_t
        # @show v_z_t
        # @show X_t

        # step
        X_vel = v_z_t * X_t - D_t * X_t
        X_tp1 = X_t + X_vel * dt_t # [g/L]
        
        sglc_vel = cglc_t * D_t - sglc_t * D_t + v_ex_glc_t * X_t
        sglc_tp1 = sglc_t + sglc_vel * dt_t  # [mM]
        
        # update
        dovPush!(SOL, "X.ts", X_tp1)            
        dovPush!(SOL, "s.EX_glc__D_e.ts", sglc_tp1)            
        dovPush!(SOL, "c.EX_glc__D_e.ts", cglc_t)   # constant
        dovPush!(SOL, "D.ts", D_t)                  # constant
        dt_tp1 = clamp(dt_max * sglc_tp1 / cglc_t, dt_min, dt_max)  # [h]
        dovPush!(SOL, "dt.ts", dt_tp1)

    end
end

