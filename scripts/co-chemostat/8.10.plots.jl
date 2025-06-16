## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: using
@time begin
    using MetX
    using MetXNetHub
    using CSV
    using DataFrames
    using CairoMakie
    using SparseArrays
    using Random
    using Bloberias
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("0.0.project.jl")
include("0.99.utils.jl")


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
let

    global lepmodel_id = [
        "iML1515.2strains",
        "PamModel-alter.2strains"
    ][2]
    global objfun_id = [
        "ΔIle.egoist.max.yield",
        "ΔMet.egoist.max.yield",
        "non-competitive.max.yield"
    ][3]

    bb = B[r"steady.state"]

    global MEANS = Dict()
    global STDS = Dict()

    for b in bb
        ctx_b = b["context", "."][]
        ctx_b["context", "objfun_id"] == objfun_id || continue
        ctx_b["context", "lepmodel_id"] == lepmodel_id || continue
        @show ctx_b["context", "lepmodel_id"]
        @show ctx_b["context", "objfun_id"]

        for (k, val) in b["stst.sol", "."]
            dovPush!(MEANS, k, val["mean"])
            dovPush!(STDS, k, val["std"])
        end
    end
    # ctx_b = bb[1]["context", "."][]
    # # ctx = ctx_b[][["context"]]
    # ctx_b[]
end

# --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
let
    
    # filter
    biomid1 = "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
    idxs1 = findall(eachindex(MEANS[biomid1])) do idx
        m = MEANS[biomid1][idx]
        s = STDS[biomid1][idx]
        abs(s) / abs(m) < 0.05
    end
    biomid2 = "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
    idxs2 = findall(eachindex(MEANS[biomid2])) do idx
        m = MEANS[biomid2][idx]
        s = STDS[biomid2][idx]
        abs(s) / abs(m) < 0.05
    end

    f = Figure()
    ax = Axis(f[1,1]; 
        title = string(
            lepmodel_id, "\n",
            objfun_id
        ), 
        xlabel = "D [1/h]",
        ylabel = "value"
    )

    yid1 = "ΔMet:X.stst"
    # yid1 = "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
    xs = MEANS["D.stst"][idxs1]
    ys = MEANS[yid1][idxs1]
    yerr = STDS[yid1][idxs1]
    scatter!(ax, xs, ys;
        markersize = 15,
        label = yid1
    )
    errorbars!(ax, xs, ys, yerr)

    yid2 = "ΔIle:X.stst"
    # yid2 = "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
    xs = MEANS["D.stst"][idxs2]
    ys = MEANS[yid2][idxs2]
    yerr = STDS[yid2][idxs2]
    scatter!(ax, xs, ys; 
        markersize = 15,
        label = yid2
    )
    errorbars!(ax, xs, ys, yerr)
    axislegend(ax)
    f
end