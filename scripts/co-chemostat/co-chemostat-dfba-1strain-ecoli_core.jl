# MARK: using
@time begin
    using EColi_Coexistence_2024
    using MetX
    using MetX.MetXNetHub
    using MetX.GLPK
    using MetX.Ipopt
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

# -- .. - .-- .-. . .... -- -- -- .. ...
#=
DOING
- Reproduce @joyStudyGrowthEscherichia2010 results
- dfba simulations of a chemostat for a two strain system.
- testing different objective functions.
- testing different initial conditions.
=# 

# -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Prepare opm
include("1.2.prepare.core_ecoli.opm.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
let
    json_file = joinpath(PaperSON_dir(), "raw.json")
    global rawData = JSON.parsefile(json_file)
    global fig1Data = rawData["data"]["Fig 1"]
    nothing
end



## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: DFBA
# dfba
SOL = Dict{String, Vector{Float64}}()
let
    # objective function
    set_linear_obj!(opm, 
        ["BIOMASS_Ecoli_core_w_GAM", "EX_glc__D_e", "EX_nh4_e", "EX_o2_e", "EX_co2_e"], 
        [1e5, -1e-2, -1e-2, -1e-2, -1e-2]
    )
    
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


## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: time plot
let
    f = Figure()
    ax = Axis(f[3:5, 1:3];
        xlabel = "time [h]"
    )

    # y_ts_id0 = "X.ts"
    # y_ts_id0 = "max_u.EX_glc__D_e.ts"
    # y_ts_id0 = "dt.ts"
    # y_ts_id0 = "s.EX_glc__D_e.ts"
    y_ts_id0 = "v.BIOMASS_Ecoli_core_w_GAM.ts"
    # y_ts_id0 = "D.ts"

    nsamples = 10_000
    idx0 = 0.004
    idx1 = 0.9

    xs = cumsum(SOL["dt.ts"])
    ys = SOL[y_ts_id0]
    idxs = _indexes(xs, ys; nsamples, idx0, idx1)
    scatter!(ax, 
        xs[idxs], ys[idxs];
        label = y_ts_id0, 
    )
    # ax.limits = (nothing, nothing, 0.0, 0.25)

    ts_id1 = "D.ts"
    xs = cumsum(SOL["dt.ts"])
    ys = SOL[ts_id1]
    idxs = _indexes(xs, ys; nsamples, idx0, idx1)
    scatter!(ax, 
        xs[idxs], ys[idxs];
        label = ts_id1
    )

    Legend(f[1:2,1:3], ax)
    f
end