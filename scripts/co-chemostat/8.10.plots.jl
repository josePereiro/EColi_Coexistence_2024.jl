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
# MARK: Colllect
let

    bb = B[r"steady.state"]

    global MEANS0 = Dict()
    global STDS0 = Dict()
    global SUMS0 = Dict()

    for (lepmodel_id, objfun_id) in Iterators.product(
        [
            "iML1515.2strains",
            "PamModel-alter.2strains"
        ],
        [
            "ΔIle.egoist.max.yield",
            "ΔMet.egoist.max.yield",
            "non-competitive.max.yield"
        ]
    )

        for b in bb
            ctx_b = b["context", "."][]
            ctx_b["context", "objfun_id"] == objfun_id || continue
            ctx_b["context", "lepmodel_id"] == lepmodel_id || continue
            @show ctx_b["context", "lepmodel_id"]
            @show ctx_b["context", "objfun_id"]
            for (k, val) in b["stst.sol", "."]
                MEANS1 = get!(MEANS0, (lepmodel_id, objfun_id), Dict())
                dovPush!(MEANS1, k, val["mean.last"])
                STDS1 = get!(STDS0, (lepmodel_id, objfun_id), Dict())
                dovPush!(STDS1, k, val["std.last"])
                SUMS1 = get!(SUMS0, (lepmodel_id, objfun_id), Dict())
                dovPush!(SUMS1, k, val["sum.tot"])
            end
        end
    end
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: GLOBALS

# MARK: ...IDS
global LEPMODEL_IDS = [
    "iML1515.2strains",
    "PamModel-alter.2strains"
]
global OBJFUN_IDS = [
    "ΔIle.egoist.max.yield",
    "ΔMet.egoist.max.yield",
    "non-competitive.max.yield"
]

# MARK: ...STST_FILTER
function STST_FILTER(MEANS1, STDS1)
    
    # filter
    I1 = eachindex(MEANS1["dt.stst"])
    return findall(I1) do idx

        m = MEANS1["ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"][idx]
        s = STDS1["ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"][idx]
        (abs(s) / abs(m) < 0.05) || return false
        
        m = MEANS1["ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"][idx]
        s = STDS1["ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"][idx]
        (abs(s) / abs(m) < 0.05) || return false

        X = MEANS1["ΔMet:X.stst"][idx]
        X > 0.1 || return false
        dt = MEANS1["dt.stst"][idx]
        dt < 0.005 || return false
        return true
    end
end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: X plots
let
    
    global lepmodel_id = [
        "iML1515.2strains",
        "PamModel-alter.2strains"
    ][1]
    global objfun_id = [
        "ΔIle.egoist.max.yield",
        "ΔMet.egoist.max.yield",
        "non-competitive.max.yield"
    ][3]

    MEANS1 = MEANS0[(lepmodel_id, objfun_id)]
    STDS1 = STDS0[(lepmodel_id, objfun_id)]

    # filter
    idxs = STST_FILTER(MEANS1, STDS1)

    f = Figure()
    ax = Axis(f[1,1]; 
        title = string(
            lepmodel_id, "\n",
            objfun_id
        ), 
        xlabel = "D [1/h]",
        ylabel = "X [g/L]", 
        limits = (0.05, 0.7, -0.05, 1.4)
    )

    yid1 = "ΔMet:X.stst"
    # yid1 = "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
    xs = MEANS1["D.stst"][idxs]
    ys = MEANS1[yid1][idxs]
    yerr = STDS1[yid1][idxs]
    scatter!(ax, xs, ys;
        markersize = 15,
        label = yid1,
        color = :red
    )
    errorbars!(ax, xs, ys, yerr; color = :red)

    yid2 = "ΔIle:X.stst"
    # yid2 = "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
    xs = MEANS1["D.stst"][idxs]
    ys = MEANS1[yid2][idxs]
    yerr = STDS1[yid2][idxs]
    scatter!(ax, xs, ys; 
        markersize = 15,
        label = yid2,
        color = :blue
    )
    errorbars!(ax, xs, ys, yerr; color = :red)

    # Add the total
    xs = MEANS1["D.stst"][idxs]
    ys = MEANS1[yid2][idxs] + MEANS1[yid1][idxs]
    yerr = STDS1[yid2][idxs] + STDS1[yid1][idxs]
    scatter!(ax, xs, ys; 
        markersize = 15,
        label = "total",
        color = :black
    )
    errorbars!(ax, xs, ys, yerr; color = :black)
    axislegend(ax)
    f
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: X1/X2 plots
let

    f = Figure()
    ax = Axis(f[1,1]; 
        title = string(
            # lepmodel_id, "\n",
            # objfun_id
        ), 
        xlabel = "D [1/h]",
        ylabel = "ΔMet:X / ΔIle:X",
        # limits = (-0.05, 0.7, -0.05, 2.0)
    )

    for (lepmodel_id, objfun_id) in Iterators.product(
        [
            "iML1515.2strains",
            "PamModel-alter.2strains"
        ],
        [
            # "ΔIle.egoist.max.yield",
            # "ΔMet.egoist.max.yield",
            "non-competitive.max.yield"
        ]
    )
        MEANS1 = MEANS0[(lepmodel_id, objfun_id)]
        STDS1 = STDS0[(lepmodel_id, objfun_id)]
        # filter
        idxs = STST_FILTER(MEANS1, STDS1)

        yid1 = "ΔMet:X.stst"
        # yid1 = "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
        xs = MEANS1["D.stst"][idxs]
        ys1 = MEANS1[yid1][idxs]
        yerr1 = STDS1[yid1][idxs]
        
        yid2 = "ΔIle:X.stst"
        ys2 = MEANS1[yid2][idxs]
        yerr2 = STDS1[yid2][idxs]
        scatter!(ax, xs, ys1 ./ ys2;
            markersize = 15,
            # label = yid1,
            # color = :red
        )
        # errorbars!(ax, xs, ys, yerr; color = :red)
        # axislegend(ax)
    end
    f
end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: sglc plots
let
    global objfun_id = [
        "non-competitive.max.yield",
        "ΔIle.egoist.max.yield",
        "ΔMet.egoist.max.yield",
    ][3]

    yid1 = "s.glc__D_e.stst"

    f = Figure()
    ax = Axis(f[1,1]; 
        title = string(
            yid1, "\n",
            objfun_id
        ), 
        xlabel = "D [1/h]",
        ylabel = "s_glc [mM]", 
        limits = (0.05, 0.7, 0.0, 4.0)
    )

    for lepmodel_id in [
            "PamModel-alter.2strains",
            "iML1515.2strains"
        ]
        # global lepmodel_id = 
        MEANS1 = MEANS0[(lepmodel_id, objfun_id)]
        STDS1 = STDS0[(lepmodel_id, objfun_id)]

        # filter
        idxs = STST_FILTER(MEANS1, STDS1)

        # yid1 = "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
        xs = MEANS1["D.stst"][idxs]
        ys = MEANS1[yid1][idxs]
        yerr = STDS1[yid1][idxs]
        scatter!(ax, xs, ys;
            markersize = 15,
            label = lepmodel_id
        )
        errorbars!(ax, xs, ys, yerr)

    end

    axislegend(ax; position = :lt)
    f
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: dt plots
let
    global objfun_id = [
        "non-competitive.max.yield",
        "ΔIle.egoist.max.yield",
        "ΔMet.egoist.max.yield",
    ][3]

    yid1 = "dt.stst"

    f = Figure()
    ax = Axis(f[1,1]; 
        title = string(
            yid1, "\n",
            objfun_id
        ), 
        xlabel = "D [1/h]",
        ylabel = "dt", 
        limits = (nothing, nothing, 0.0, 60.0)
    )

    for lepmodel_id in [
            "PamModel-alter.2strains",
            "iML1515.2strains"
        ]
        # global lepmodel_id = 
        MEANS1 = MEANS0[(lepmodel_id, objfun_id)]
        STDS1 = STDS0[(lepmodel_id, objfun_id)]
        SUMS1 = SUMS0[(lepmodel_id, objfun_id)]

        # filter
        idxs = STST_FILTER(MEANS1, STDS1)

        # yid1 = "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
        xs = MEANS1["D.stst"][idxs]
        ys = SUMS1[yid1][idxs]
        scatter!(ax, xs, ys;
            markersize = 15,
            label = lepmodel_id
        )
        # errorbars!(ax, xs, ys, yerr)

    end

    # axislegend(ax; position = :lt)
    f
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: tot yield plot
let
    global objfun_id = [
        "non-competitive.max.yield",
        "ΔIle.egoist.max.yield",
        "ΔMet.egoist.max.yield",
    ][3]

    f = Figure()
    ax = Axis(f[1,1]; 
        title = string(
            "X.glc.yield", "\n",
            objfun_id
        ), 
        xlabel = "D [1/h]",
        ylabel = "Y_Xglc [g/mmol]", 
        limits = (0.05, 0.7, 0.08, 0.13)
    )

    for lepmodel_id in [
            "PamModel-alter.2strains",
            "iML1515.2strains"
        ]
        # global lepmodel_id = 
        MEANS1 = MEANS0[(lepmodel_id, objfun_id)]
        STDS1 = STDS0[(lepmodel_id, objfun_id)]

        # filter
        idxs = STST_FILTER(MEANS1, STDS1)

        # TODO: Ask chatGPT
        # yield = \mu / u_glc
        # X * D / (c - s) * D
        # [g/Lh] / (mmol/Lh)
        # [g/mmol]
        xs = MEANS1["D.stst"][idxs]
        totX = MEANS1["ΔMet:X.stst"][idxs] .+ MEANS1["ΔIle:X.stst"][idxs]
        dglc = (MEANS1["c.glc__D_e.stst"][idxs] .- MEANS1["s.glc__D_e.stst"][idxs]) 
        ys = totX ./ dglc
        # yerr = STDS1[yid1][idxs]
        scatter!(ax, xs, ys;
            markersize = 15,
            label = lepmodel_id
        )
        # errorbars!(ax, xs, ys, yerr)

    end

    axislegend(ax; position = :rt)
    f
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: strain yield plot
let
    global objfun_id = [
        "non-competitive.max.yield",
        "ΔIle.egoist.max.yield",
        "ΔMet.egoist.max.yield",
    ][3]

    global strainID = [
        "ΔIle", "ΔMet"
    ][2]

    f = Figure()
    ax = Axis(f[1,1]; 
        title = string(
            "X.glc.yield", "\n",
            objfun_id
        ), 
        xlabel = "D [1/h]",
        ylabel = "Y_Xglc [g/mmol]", 
        limits = (0.05, 0.7, 0.08, 0.13)
    )

    for lepmodel_id in [
            "PamModel-alter.2strains",
            "iML1515.2strains"
        ]
        # global lepmodel_id = 
        MEANS1 = MEANS0[(lepmodel_id, objfun_id)]
        STDS1 = STDS0[(lepmodel_id, objfun_id)]

        # filter
        idxs = STST_FILTER(MEANS1, STDS1)

        # TODO: Ask chatGPT
        # yield = D / u_glc
        # [g\Lh] / [mmol/ gh]
        # [g\L] / [mmol/ g]
        xs = MEANS1["D.stst"][idxs]
        totX = MEANS1["$strainID:BIOMASS_Ec_iML1515_WT_75p37M_f"][idxs]
        dglc = MEANS1["$strainID:v:EX_glc__D_e_b.stst"][idxs] 
        ys = totX ./ abs.(dglc)
        # yerr = STDS1[yid1][idxs]
        scatter!(ax, xs, ys;
            markersize = 15,
            label = string("$strainID:", lepmodel_id)
        )

        # xs = MEANS1["D.stst"][idxs]
        # totX = MEANS1["ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"][idxs]
        # dglc = MEANS1["ΔMet:v:EX_glc__D_e_b.stst"][idxs] 
        # ys = totX ./ abs.(dglc)
        # # yerr = STDS1[yid1][idxs]
        # scatter!(ax, xs, ys;
        #     markersize = 15,
        #     label = string("ΔMet:", lepmodel_id),
        #     marker = :rect
        # )
        # errorbars!(ax, xs, ys, yerr)

    end

    axislegend(ax; position = :rt)
    f
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: glucose intake
let
    global objfun_id = [
        "non-competitive.max.yield",
        "ΔIle.egoist.max.yield",
        "ΔMet.egoist.max.yield",
    ][3]

    f = Figure()
    ax = Axis(f[1,1]; 
        title = string(
            "X.glc.yield", "\n",
            objfun_id
        ), 
        # xlabel = "D [1/h]",
        # ylabel = "Y_Xglc [g/mmol]", 
        # limits = (0.05, 0.7, 0.08, 0.13)
    )

    for lepmodel_id in [
            "PamModel-alter.2strains",
            "iML1515.2strains"
        ]
        # global lepmodel_id = 
        MEANS1 = MEANS0[(lepmodel_id, objfun_id)]
        STDS1 = STDS0[(lepmodel_id, objfun_id)]

        # filter
        idxs = STST_FILTER(MEANS1, STDS1)

        # (c - s) * D / X
        xs = MEANS1["ΔIle:v:EX_glc__D_e_b.stst"][idxs] .* MEANS1["ΔIle:X.stst"][idxs]
        ys = (MEANS1["c.glc__D_e.stst"][idxs] .- MEANS1["s.glc__D_e.stst"][idxs]) .* MEANS1["D.stst"][idxs] .* MEANS1["ΔIle:X.stst"][idxs] ./ (MEANS1["ΔIle:X.stst"][idxs] .+ MEANS1["ΔMet:X.stst"][idxs] )
        scatter!(ax, abs.(xs), abs.(ys);
            markersize = 15,
            label = lepmodel_id
        )

    end

    axislegend(ax; position = :lt)
    f
end