# MARK: using
@time begin
    using EColi_Coexistence_2024
    using MetX
    using MetXNetHub
    using GLPK
    using Ipopt
    using Clp
    using CSV
    using DataFrames
    using JSON
    using CairoMakie
    using Serialization
    using SparseArrays
    using Random
    using Statistics
end

# -- .. - .-- .-. . .... -- -- -- .. ...
include("0.0.project.jl")
include("0.99.utils.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
# LOAD
BLOBS = []
DATDIR="/Users/Pereiro/.julia/dev/EColi_Coexistence_2024/scripts/co-chemostat/data/3.1.co-chemostat-dfba-2strain-PAM_model"
for file in readdir(DATDIR; join = true)
    endswith(file, ".jls") || continue
    @show file
    SOL = deserialize(file)
    push!(BLOBS, SOL)
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Exploring
let
    for blob in BLOBS
        println("--------")
        @show blob["scope"][:dt_t0]
        @show blob["scope"][:D_t0]
        @show blob["scope"][:prepare_lepmodel]
        @show blob["scope"][:objfun_id]
    end
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Utils
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
    if nsamples > 0 && nidxs > nsamples
        return rand(idxs, nsamples)
    end
    return idxs
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: X vs D
let
    objfun_id = "ΔIle.egoist.max.yield"
    # objfun_id = "ΔMet.egoist.max.yield"
    # prepare_lepmodel = "iML1515.2strains"
    prepare_lepmodel = "PamModel-alter.2strains"

    xs = Float64[]
    ΔIle_ys = Float64[]
    ΔMet_ys = Float64[]
    for blob in BLOBS
        scope = blob["scope"]
        sol = blob["SOL"]

        # FILTER
        scope[:objfun_id] == objfun_id || continue
        scope[:prepare_lepmodel] == prepare_lepmodel || continue
        
        ΔIle_y = mean(last(sol["ΔIle:X.ts"], 10))
        ΔMet_y = mean(last(sol["ΔMet:X.ts"], 10))

        push!(xs, scope[:D_t0])
        push!(ΔMet_ys, ΔMet_y)
        push!(ΔIle_ys, ΔIle_y)
    end

    f = Figure()

    ax = Axis(f[3:5, 1:3];
        title = string(objfun_id, "\n", prepare_lepmodel),
        xlabel = "D [h]",
        ylabel = "X [g/L]",
        limits = (nothing, nothing, 0.0, 1.0)
    )

    scatter!(ax, 
        xs, ΔIle_ys;
        label = "ΔIle:X", 
    )
    scatter!(ax, 
        xs, ΔMet_ys;
        label = "ΔMet:X", 
    )
    Legend(f[1:2,1:3], ax)
    f
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: dynamics
let
    PKG = BLOBS[1]
    scope = PKG["scope"]
    SOL = PKG["SOL"]

    f = Figure()
    ax = Axis(
        f[3:5, 1:3];
        title = string(
            "D=", scope[:D_t0], "\n", 
            "dt_t0: ", scope[:dt_t0], "\n",
            scope[:objfun_id], "\n",
            scope[:prepare_lepmodel],
        ),
        xlabel = "time [h]",
        # limits = (nothing, nothing, 0.0, 1.0)
    )

    # y_ts_id0 = "ΔMet:v:EX_glc__D_e_b.ts"
    # y_ts_id1 = "ΔIle:v:EX_glc__D_e_b.ts"
    # y_ts_id0 = "ΔMet:X.ts"
    # y_ts_id1 = "ΔIle:X.ts"
    # y_ts_id0 = "max_u.glc__D_e.ts"
    # y_ts_id0 = "dt.ts"
    # y_ts_id0 = "s.glc__D_e.ts"
    # "EX_ile__L_e", "EX_met__L_e"
    y_ts_id0 = "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
    y_ts_id1 = "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
    # y_ts_id0 = "D.ts"
    # y_ts_id1 = ""
    @show y_ts_id0
    @show y_ts_id1

    nsamples = 5_000
    idx0 = 0.004
    idx1 = 0.9

    for datid in [y_ts_id0, y_ts_id1]
        isempty(datid) && continue
        isnothing(datid) && continue
        
        xs = cumsum(SOL["dt.ts"])
        ys = SOL[datid]
        idxs = _indexes(xs, ys; nsamples, idx0, idx1)
        scatter!(ax, 
            xs[idxs], ys[idxs];
            label = datid, 
        )
    end

    Legend(f[1:2,1:3], ax)
    f
end