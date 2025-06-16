# MARK: using
@time begin
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

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("0.0.project.jl")
include("0.99.utils.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
let
    bb = blobbatch(B, "co-chemostat-fba-2strains-iML1515-based")

    # lepmodel_id = "iML1515.2strains"
    # objfun_id = "non-competitive.max.yield"

    rxns_rawids =  ["EX_glc__D_e", "EX_met__L_e", "EX_ile__L_e", "BIOMASS_Ec_iML1515_WT_75p37M", "EX_o2_e", "EX_co2_e"]
    
    global DAT = Dict()
    for b in bb
        lepmodel_id = b["context", "lepmodel_id"]
        objfun_id = b["context", "objfun_id"]
        sol = b["sol", "val"]
        D = b["context", "D"]

        dat = get!(DAT, (objfun_id, lepmodel_id), Dict())
        dovPush!(dat, "D", D)

        for rawid in rxns_rawids
            # ΔIle
            flux_b = sol["ΔIle:$(rawid)_b"]
            flux_f = sol["ΔIle:$(rawid)_f"]
            dovPush!(dat, "ΔIle:$(rawid)", abs(flux_f - flux_b))
            # ΔMet
            flux_b = sol["ΔMet:$(rawid)_b"]
            flux_f = sol["ΔMet:$(rawid)_f"]
            dovPush!(dat, "ΔMet:$(rawid)", abs(flux_f - flux_b))
            # Medium
            if (rawid != "BIOMASS_Ec_iML1515_WT_75p37M")
                flux_b = sol["Medium:$(rawid)_b"]
                flux_f = sol["Medium:$(rawid)_f"]
                dovPush!(dat, "Medium:$(rawid)", abs(flux_f - flux_b))
            end
        end
    end
    DAT
end


## -- .. - .-- .-. . .... -- -- -- .. ...
function _normalizare!(vec) 
    vec .= vec ./ maximum(vec)
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Plot 1
let
    rxns_rawids =  ["EX_glc__D_e", "EX_met__L_e", "EX_ile__L_e", "BIOMASS_Ec_iML1515_WT_75p37M", "EX_o2_e", "EX_co2_e"]

    lepmodel_id = "iML1515.2strains"
    
    objfun_ids = [
        "ΔIle.egoist.max.yield",
        "ΔMet.egoist.max.yield",
        "ΔIle.egoist.only",
        "ΔMet.egoist.only",
        "non-competitive.max.yield" 
    ]
    dat = DAT[(objfun_id, lepmodel_id)]
    
    raw1 = rxns_rawids[6]
    raw2 = rxns_rawids[5]
    
    f = Figure(; size = (600,600))
    ax = Axis(f[3:5,1:3];
        title = join([objfun_id, lepmodel_id], "\n"),
        xlabel = "D",
        ylabel = "value",
    )

    xs = dat["D"]
    yid = "ΔMet:$(raw1)"
    ys = dat[yid]
    scatter!(xs, ys; label = yid)
    
    xs = dat["D"]
    yid = "ΔIle:$(raw2)"
    ys = dat[yid]
    scatter!(xs, ys; label = yid)

    Legend(f[1:2,1:3], ax)
    f
end


## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Plot 1
let
    f = Figure(; size = (600,600))
    
    ΔMetRxn = "EX_glc__D_e"
    # ΔMetRxn = "EX_met__L_e"
    ΔIleRxn = "EX_glc__D_e"
    # ΔIleRxn = "EX_ile__L_e"
    
    ax = Axis(f[3:5,1:3];
        xlabel = "D [1/h]",
        ylabel = "value",
        aspect=1.5
    )
    for (simid, color) in [
        ("ΔMet:max.glc.yield", :blue),
        ("ΔIle:max.glc.yield", :red)
    ]
        ys = 
            SOL["$simid:D"] ./
            abs.(
                SOL["$simid:ΔMet:$(ΔMetRxn)_f"] - 
                SOL["$simid:ΔMet:$(ΔMetRxn)_b"]
            )
        scatter!(ax, SOL["$simid:D"],
            # _normalizare!(ys);
            ys;
            label = "$simid - D/ΔMet:$(ΔMetRxn)",
            marker = :circle,
            color
        )
        ys = 
            SOL["$simid:D"] ./
            abs.(
                SOL["$simid:ΔIle:$(ΔIleRxn)_f"] - 
                SOL["$simid:ΔIle:$(ΔIleRxn)_b"]
            )
        scatter!(ax, SOL["$simid:D"],
            # _normalizare!(ys);
            ys;
            label = "$simid - D/ΔIle:$(ΔIleRxn)",
            marker = :rect,
            color
        )
    end
    Legend(f[1:2,1:3], ax)
    f
end
