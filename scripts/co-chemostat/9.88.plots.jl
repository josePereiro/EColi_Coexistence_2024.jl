## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
@time begin
    using CSV
    using MetX
    using MetX: Clp, GLPK, Ipopt
    using JSON
    using Random
    using MetXNetHub
    using DataFrames
    using CairoMakie
    using SparseArrays
    using Serialization
    using EColi_Coexistence_2024
    using Scoperias
    using Bloberias
    using Dates
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("0.0.project.jl")
include("0.99.utils.jl")
include("9.99.utils.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# TODO: create SimuleOs


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: biom heatmap
let
    bb =  B[r"fba.glc.q_aa.phy"]

    # filter control
    aux_id = "ΔIle"

    xs, ys, zs = Float64[], Float64[], Float64[]
    
    for b in bb

        b["context", "aux_id"] == aux_id || continue

        aa_U = b["context", "aa_U"]
        glc_U = b["context", "glc_U"]
        aux_ex = b["context", "aux_ex"]
        glc_ex = b["context", "glc_ex"]

        opm_sol = b["opm.sol", "."]
        biom = opm_sol["BIOMASS_Ec_iML1515_WT_75p37M_f"]
        glc_q = opm_sol[glc_ex]
        aa_q = opm_sol[aux_ex]

        push!(xs, glc_U)
        push!(ys, aa_U)
        push!(zs, biom)
        
    end

    # Create figure and axis
    f = Figure()
    ax = Axis(f[1, 1];
        title = aux_id,
        xlabel = "U_glc (mmol/gh)",
        ylabel = "U_aa (mmol/gh)",
        aspect = 1.0
    )

    # Plot scatter with color mapped to zs
    sc = scatter!(ax, xs, ys; 
        color=zs,
        marker = :rect,
        markersize = 27
    )

    # Add colorbar, linked to the scatter object
    Colorbar(f[1, 2], sc;
        label = "growth rate (1/h)"
    )

    f
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: q lines
let
    bb =  B[r"fba.glc.aa.phy"]

    # filter control
    aux_id = "ΔIle"

    U_LINES = U_lines(aux_id)

    nonlim_arg = "aa_U"

    if nonlim_arg == "glc_U"
        xname = "aa_U"
    elseif nonlim_arg == "aa_U"
        xname = "glc_U"
    end

    f = Figure()
    ax = Axis(f[1,1];
        title = aux_id,
        xlabel = "$(xname) (mmol/gh)",
        ylabel = "growth rate (1/h)",
        limits = (nothing, nothing, nothing, 0.8)
    )

    Umax = maximum(keys(U_LINES[nonlim_arg]))
    for (U0, line) in U_LINES[nonlim_arg]
        idxs = sortperm(line["xs"])
        lines!(ax, line["xs"][idxs], line["ys"][idxs];
            color = :black,
            alpha = 0.7, 
            linewidth = max(1, 5 * U0 / Umax)
        )
    end

    # Env
    ENV = U_envelop(U_LINES, nonlim_arg)

    scatter!(ax, ENV[:xs], ENV[:ys];
        marker = :circle, 
        markersize = 20, 
        color = :black
    )

    # write envelop

    f
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
let
    a = 1
    nm = dfname("growth", (;a))
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: s lines
let
    bb =  B[r"fba.glc.aa.phy"]

    # filter control
    aux_id = "ΔIle"

    KM_glc = 1.75 # mM
    Vmax_glc = 10.0   # mmol/gh

    if aux_id == "ΔMet"
        Vmax_aa = 0.4   # mmol/gh
        KM_aa = 0.2     # mM
    elseif aux_id == "ΔIle"  
        Vmax_aa = 0.4   # mmol/gh
        KM_aa = 0.2     # mM
    end

    U_LINES = U_lines(aux_id)
    
    nonlim_arg = "aa_U"

    if nonlim_arg == "glc_U"
        xname = "aa_s"
    elseif nonlim_arg == "aa_U"
        xname = "glc_s"
    end

    f = Figure()
    ax = Axis(f[1,1];
        title = aux_id,
        xlabel = "$(xname) (mM)",
        ylabel = "growth rate (1/h)",
        limits = (nothing, nothing, nothing, 0.8)

    )

    glc_U1 = maximum(keys(U_LINES[nonlim_arg]))
    for (glc_U, line) in U_LINES[nonlim_arg]
        idxs = sortperm(line["xs"])
        lines!(ax, 
            line["xs"][idxs], line["ys"][idxs];
            color = :black,
            alpha = 0.7, 
            linewidth = max(1, 5 * glc_U / glc_U1)
        )
    end

    # Env
    ENV = U_envelop(U_LINES, nonlim_arg)
    scatter!(ax, ENV[:xs], ENV[:ys];
        marker = :circle, 
        markersize = 20, 
        color = :black
    )

    f

end


