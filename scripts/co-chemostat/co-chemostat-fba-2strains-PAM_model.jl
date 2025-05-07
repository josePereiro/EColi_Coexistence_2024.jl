# MARK: using
@time begin
    using MetX
    using MetX.MetXNetHub
    using MetX.GLPK
    using MetX.Ipopt
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
Simulate a chemostat for two strains at different dilution rates
    See if it is the case that the dilution rate is not affecting coexistance
=# 

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Prepare opm
include("1.1.prepare.PamModel-alter.2strains.opm.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Compute D serie
targetDs = shuffle(collect(range(0.1, 0.7; length=50)))
# set_linear_obj!(opm, 
#     ["ΔMet:EX_glc__D_e_b", "ΔIle:EX_glc__D_e_b"], 
#     [-1, -1]
# )
SOL = Dict()
for (simid, D) in Iterators.product(
        [
            "ΔMet:max.glc.yield",
            "ΔIle:max.glc.yield",
        ],
        targetDs
    )

    # SETUP
    if (simid == "ΔMet:max.glc.yield")
        c = [-1e3, -1, -1, -1, -1, -1, -1]
    elseif (simid == "ΔIle:max.glc.yield")
        c = [-1, -1e3, -1, -1, -1, -1, -1]
    end
    set_linear_obj!(opm, [
            "ΔMet:EX_glc__D_e_b", 
            "ΔIle:EX_glc__D_e_b",
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ], 
        c
    )

    _eps = 1e-3
    bounds!(opm, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f", D - _eps, D + _eps)
    bounds!(opm, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f", D - _eps, D + _eps)

    # OPTIMIZED
    try
        @time optimize!(opm)
    catch e
        (e isa InterruptException) && break
        continue
    end

    # COLLECT
    dovPush!(SOL, "$simid:D", D)
    for (comp, flux) in Iterators.product(
            ["ΔMet", "ΔIle", "Medium"],
            [
                "EX_glc__D_e_f",
                "EX_glc__D_e_b",
                "EX_met__L_e_f",
                "EX_met__L_e_b",
                "EX_ile__L_e_f",
                "EX_ile__L_e_b",
                "BIOMASS_Ec_iML1515_WT_75p37M_f",
                "BIOMASS_Ec_iML1515_WT_75p37M_b"
            ]
        )

        try
            dovPush!(SOL, "$simid:$comp:$flux", 
                solution(opm, "$comp:$flux")
            )
        catch; end
    end

    # INFO
    println("="^30)
    @show D
    @show simid
    @show objective_value(opm)
    
    _opm_summary_PamModel_alter(opm)

    println("-"^30)
    @show D / solution(opm, "ΔMet:EX_glc__D_e_b") 
    @show D / solution(opm, "ΔIle:EX_glc__D_e_b") 
    
end

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Plots

## -- .. - .-- .-. . .... -- -- -- .. ...
function _normalizare!(vec) 
    vec .= vec ./ maximum(vec)
end

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

## -- .. - .-- .-. . .... -- -- -- .. ...
let
end

## -- .. - .-- .-. . .... -- -- -- .. ...
