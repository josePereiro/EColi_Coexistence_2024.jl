# MARK: using
@time begin
    using EColi_Coexistence_2024
    using MetX
    using MetXNetHub
    using GLPK
    using Ipopt
    using CSV
    using DataFrames
    using JSON
    using CairoMakie
    using Serialization
    using SparseArrays
    using Random
end

## -- .. - .-- .-. . .... -- -- -- .. ...
let
    json_file = joinpath(PaperSON_dir(), "raw.json")
    rawData = JSON.parsefile(json_file)
    fig1Data = rawData["data"]["Fig 1"]
    nothing
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