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
# DOING detect stst
let
    # Random.seed!(123)
    global _b = rand(rand(B, r"co-"));
    _v = _b["dyn.sol", "."]["ΔIle:X.ts"];
    w = 0.1
    i0, i1 = firstindex(_v), lastindex(_v)
    l = length(_v)
    i = ceil(Int, l * (1 - w))
    i = clamp(i, i0, i1)
    m, s = mean(_v[i:i1]), std(_v[i:i1])
    @show m, s
    s / abs(m) < 0.1
end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# TODO/ Move to Bloberias

function stateful_blobs(B::Bloberia, bbid_prefix = nothing)
    return Channel{bBlob}(1) do _ch
        bbs = eachbatch(B, bbid_prefix;
            sortfun = sort!,
            ch_size = 1,
            n_tasks = 1
        )
        for bb in bbs
            for b in bb
                put!(_ch, b)
            end
        end
    end
end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
_blob_CH = stateful_blobs(B,
    r"co-chemostat-dfba-2strains-iML1515-based"
)

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
let
    # dyn.sol
    #  "ΔIle:X.ts"
    #   "dt.ts"
    #   "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
    #   "ΔMet:U:EX_glc__D_e_b.ts"
    #   "c.glc__D_e.ts"
    #   "ΔIle:v:EX_glc__D_e_b.ts"
    #   "ΔIle:U:EX_glc__D_e_b.ts"
    #   "ΔMet:X.ts"
    #   "ΔMet:v:EX_glc__D_e_b.ts"
    #   "s.glc__D_e.ts"
    #   "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
    #   "D.ts"

    # context
    # ["st_label_0x94cca352d2e8e22d", "C", "blockuuid", "c_glc__D_e_t0", "litecontext", "lepmodel_id", "Main", "prepare_iML1515_2strains_opm", "D_t0", "prepare_PamModel_alter_2strains_opm", "B", "blockbatch!", "prepare_lepmodel_ref", "select_objfun", "SOL", "dovLast", "CLI", "double_model", "locate", "net0", "glc_Km", "glc_Vmax", "PROJ", "@cli", "Base", "prepare_lepmodel", "X_min", "niters", "lep0", "@init_block", "blockttag", "@blockbatch!", "ascontext", "Core", "dovLast!", "dt_t0", "bb", "lepmodel_args", "s_glc__D_e_t0", "objfun_id", "dt_min", "objiders", "PROJVER", "dovPush!", "ΔMet_X_t0", "u_max", "G", "opm", "ΔIle_X_t0", "X_max", "dt_max", "st_label_0x914465398f82596b", "blockver", "prepare_ecoli_core_opm", "objcoes"]

    # bb = rand(B, r"co-chemostat-dfba-2strains-iML1515-based")
    # b = rand(bb)

    # @show bb.id
    # @show b[]
    b = take!(_blob_CH)
    
    lepmodel_id = b["context", "lepmodel_id"]
    objfun_id = b["context", "objfun_id"]
    D_t0 = b["context", "D_t0"]

    f = Figure(; size = (600,600))
    ax1 = Axis(f[3:5,1:3];
        title = join([
            objfun_id, 
            lepmodel_id, 
            "D = $(D_t0) [1/h]"
        ], "\n"),
        xlabel = "time [h]",
        ylabel = "value1",
        limits = (nothing, nothing, 0.0, nothing)
    )

    ax2 = Axis(f[3:5,1:3]; 
        ylabel="value2", 
        yaxisposition = :right, 
        xticklabelsvisible=false, 
        xgridvisible=false, 
        ygridvisible=false,
        spinewidth=0, 
    )

    dsol = b["dyn.sol", "."]
    xs = cumsum(dsol["dt.ts"])
    
    yid1 = "ΔIle:X.ts"
    ys = dsol[yid1]
    sc1 = scatter!(ax1, xs, ys; color = :blue)
    
    yid2 = "ΔMet:X.ts"
    ys = dsol[yid2]
    sc2 = scatter!(ax1, xs, ys; color = :red)

    yid3 = "s.glc__D_e.ts"
    ys = dsol[yid3]
    sc3 = scatter!(ax2, xs, ys; color = :black)

    # Legend(f[1:2,1:3], [sc1, sc2], [yid1, yid2])
    Legend(f[1:2,1:3], [sc1, sc2, sc3], [yid1, yid2, yid3])
    return f
    
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
let
    global bb = blobbatch(B, "co-chemostat-dfba-2strains-iML1515-based")
    nbs = length(bb)

    # dyn.sol
    #  "ΔIle:X.ts"
    #   "dt.ts"
    #   "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f"
    #   "ΔMet:U:EX_glc__D_e_b.ts"
    #   "c.glc__D_e.ts"
    #   "ΔIle:v:EX_glc__D_e_b.ts"
    #   "ΔIle:U:EX_glc__D_e_b.ts"
    #   "ΔMet:X.ts"
    #   "ΔMet:v:EX_glc__D_e_b.ts"
    #   "s.glc__D_e.ts"
    #   "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
    #   "D.ts"


    # context
    # ["st_label_0x94cca352d2e8e22d", "C", "blockuuid", "c_glc__D_e_t0", "litecontext", "lepmodel_id", "Main", "prepare_iML1515_2strains_opm", "D_t0", "prepare_PamModel_alter_2strains_opm", "B", "blockbatch!", "prepare_lepmodel_ref", "select_objfun", "SOL", "dovLast", "CLI", "double_model", "locate", "net0", "glc_Km", "glc_Vmax", "PROJ", "@cli", "Base", "prepare_lepmodel", "X_min", "niters", "lep0", "@init_block", "blockttag", "@blockbatch!", "ascontext", "Core", "dovLast!", "dt_t0", "bb", "lepmodel_args", "s_glc__D_e_t0", "objfun_id", "dt_min", "objiders", "PROJVER", "dovPush!", "ΔMet_X_t0", "u_max", "G", "opm", "ΔIle_X_t0", "X_max", "dt_max", "st_label_0x914465398f82596b", "blockver", "prepare_ecoli_core_opm", "objcoes"]

    
    # b = bb[rand(1:nbs)]
    b = bb[48]

    lepmodel_id = b["context", "lepmodel_id"]
    objfun_id = b["context", "objfun_id"]
    D_t0 = b["context", "D_t0"]
    # D_t0 = b["context", "dt_t0"]

    f = Figure(; size = (600,600))
    ax1 = Axis(f[3:5,1:3];
        title = join([
            objfun_id, 
            lepmodel_id, 
            "D = $(D_t0) [1/h]"
        ], "\n"),
        xlabel = "time [h]",
        ylabel = "value1",
        limits = (nothing, nothing, 0.0, 1.3)
    )

    # ax2 = Axis(f[3:5,1:3]; 
    #     yaxisposition=:right, 
    #     ylabel="value2", 
    #     xticksvisible=false, 
    #     xgridvisible=false, 
    #     ygridvisible=false
    # )

    dsol = b["dyn.sol", "."]
    xs = cumsum(dsol["dt.ts"])
    
    yid1 = "ΔIle:X.ts"
    ys = dsol[yid1]
    sc1 = scatter!(ax1, xs, ys; color = :blue)
    
    yid2 = "ΔMet:X.ts"
    ys = dsol[yid2]
    sc2 = scatter!(ax1, xs, ys; color = :red)

    yid3 = "dt.ts"
    ys = dsol[yid3]
    # sc3 = scatter!(ax2, xs, ys; color = :black)

    Legend(f[1:2,1:3], [sc1, sc2], [yid1, yid2])
    # Legend(f[1:2,1:3], [sc1, sc2, sc3], [yid1, yid2, yid3])
    f

end

## --.. . -. -. .- . . .. - - --  . .. . .- - . << >> -. -
let

    # Sample data
    x = LinRange(0, 2π, 500)
    y1 = sin.(x)
    y2 = cos.(x)

    # Create a figure
    fig = Figure()

    # Create two axes: one on top of the other
    ax1 = Axis(fig[1, 1], ylabel="sin(x)", xlabel="x")
    ax2 = Axis(fig[1, 1], yaxisposition=:right, ylabel="cos(x)", xticksvisible=false, xgridvisible=false, ygridvisible=false)

    # Plot sin(x) on the left y-axis
    lines!(ax1, x, y1, color=:blue, label="sin(x)")

    # Plot cos(x) on the right y-axis
    lines!(ax2, x, y2, color=:red, label="cos(x)")

    # Optional: add legends manually if desired
    axislegend(ax1, position=:lt)
    axislegend(ax2, position=:rt)

    fig

end