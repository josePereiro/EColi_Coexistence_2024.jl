## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: using
@time begin
    using EColi_Coexistence_2024
    using ProjFlows
    using CSV
    using JSON
    using DataFrames
    using CairoMakie
    using Serialization
    using Scoperias
    using JuMP
    using Ipopt
    using MetX
    using NDHistograms
    using NDHistograms: descretize
    using Statistics
end

## . -. - -.-- . . - .- .--. - < / -. .- - 
## MARK: get biomass composition
## 
## 

let

    global net0 = pull_net("iML1515")
    biom = extras(net0, "BIOM")
    global biom_comp = Dict{String, Float64}()
    for (i, s) in enumerate(stoi(net0, :, biom))
        s < 0 || continue
        biom_comp[net0.mets[i]] = abs(s) # mmol/gCDW
    end

    # BIOMASS
    # ["murein3p3p_p", "10fthf_c", "hemeO_c", "pe160_p", "lys__L_c", "nadh_c", "ptrc_c", "his__L_c", "ni2_c", "colipa_e", "pg161_p", "chor_c", "pheme_c", "succoa_c", "sheme_c", "ctp_c", "dttp_c", "udcpdp_c", "mocogdp_c", "zn2_c", "utp_c", "spmd_c", "dgtp_c", "gtp_c", "ile__L_c", "amet_c", "2dmmql8_c", "enter_c", "accoa_c", "gly_c", "tyr__L_c", "clpn181_p", "thr__L_c", "mg2_c", "pydx5p_c", "murein3px4p_p", "mql8_c", "ribflv_c", "dctp_c", "fad_c", "atp_c", "adocbl_c", "pro__L_c", "murein4px4px4p_p", "cu2_c", "phe__L_c", "thf_c", "murein4p4p_p", "val__L_c", "mlthf_c", "btn_c", "datp_c", "leu__L_c", "so4_c", "4fe4s_c", "glu__L_c", "mn2_c", "coa_c", "nad_c", "k_c", "pe161_p", "ser__L_c", "ca2_c", "nh4_c", "clpn161_p", "malcoa_c", "arg__L_c", "5mthf_c", "glycogen_c", "cl_c", "gln__L_c", "trp__L_c", "mobd_c", "q8h2_c", "thmpp_c", "mococdp_c", "nadph_c", "fe2_c", "ala__L_c", "cys__L_c", "2fe2s_c", "pg160_p", "fe3_c", "clpn160_p", "asn__L_c", "pg181_p", "h2o_c", "murein4px4p_p", "pe181_p", "cobalt2_c", "lipopb_c", "nadp_c", "asp__L_c", "gthrd_c", "met__L_c"]

end

## . -. - -.-- . . - .- .--. - <\ > -. .- - 
## MARK: wintermute
## 
## 

## . -. - -.-- . . - .- .--. - < /> -. .- - 
function load_sheet(datasheet)

    psondir = EColi_Coexistence_2024.PaperSON_dir()
    winterdir = joinpath(psondir, "data/wintermuteEmergentCooperationMicrobial2010")
    raw = JSON.parsefile(joinpath(winterdir, "raw.json"))
    
    supraw = JSON.parsefile(joinpath(winterdir,
        raw["data"]["msb201066-sup-0002.xls"]["PaperSON.file"]
    ))
    
    fn = joinpath(winterdir, supraw["data"][datasheet]["data"]["PaperSON.file"])
    df = CSV.read(fn, DataFrame; header = 2)
    select!(df, Not(1:2))
    return df
end

function plot_mat(df)

    # Example data
    # Short text labels
    tickslabels = names(df)
    A = Matrix(df)

    # Create figure and axis
    fig = Figure(;size = (900,900))
    ax = Axis(fig[1,1],
        xticks = (eachindex(tickslabels), tickslabels),
        yticks = (eachindex(tickslabels), tickslabels),
        xticklabelrotation = π/2,
        yticklabelrotation = 0,
    )

    # Plot heatmap
    hm = heatmap!(ax, A)

    # Add colorbar linked to heatmap
    Colorbar(fig[1,2], hm)

    # Display or save
    fig

end

# . -. - -.-- . . - .- .--. - < / -. .- - 
let
    global N1_mCherry = load_sheet("mCherry_Counts")
    global N1_Venus = load_sheet("Venus_Counts")
    N1_all = N1_mCherry .+ N1_Venus

    X0 = load_sheet("Day1_Rep1")
    X1 = load_sheet("Day4_Rep1")

    N0_all = (N1_all .* X0) ./ X1
    N0_mCherry = (N0_all ./ 2)
    N0_Venus = (N0_all ./ 2)
    
    global F1_mCherry = clamp.(N1_mCherry ./ N0_mCherry, 0, 1e7)
    global F1_Venus = clamp.(N1_Venus ./ N0_Venus, 0, 1e7)

    global ΔN_mCherry = clamp.(N1_mCherry .- N0_mCherry, 0.0, 1e7)
    global ΔN_Venus = clamp.(N1_Venus .- N0_Venus, 0.0, 1e7)

    # plot_mat(clamp.(ΔN_mCherry .+ ΔN_Venus, 0.0, 1e7))
    plot_mat(F1_mCherry .+ F1_Venus)
end

## . -. - -.-- . . - .- .--. - < / -. .- - 
# MARK: MUTANT MAP
let
    global AUX_MAP = Dict(
        "argB" => "arg__L_c",
        "argC" => "arg__L_c",
        "argE" => "arg__L_c",
        "argH" => "arg__L_c",
        "proC" => "pro__L_c",
        "icd" => "glu__L_c",
        "nadC" => "nac_c",      # "nicotinate (niacin) / NAD precursor",
        "panB" => "pnto__R_c",  # "pantothenate (vitamin B5)",
        "panD" => "pnto__R_c",  # "pantothenate (vitamin B5)",
        "serC" => "ser__L_c",
        "cysC" => "cys__L_c", # "cysteine (or a sulfide source)",
        "cysE" => "cys__L_c", # "cysteine",
        "cysH" => "cys__L_c", # "cysteine or possibly sulfite",
        "glnA" => "gln__L_c", # "glutamine",
        "serA" => "ser__L_c", # "serine",
        "hisB" => "his__L_c", # "histidine",
        "hisC" => "his__L_c", # "histidine",
        "hisD" => "his__L_c", # "histidine",
        "hisI" => "his__L_c", # "histidine",
        "metA" => "met__L_c", # "methionine",
        "guaA" => "gua_c", # "guanine (GMP precursor)",
        "guaB" => "gua_c", # "guanine (XMP/ GMP precursor)",
        "purA" => "ade_c", # "adenine (AMP precursor)",
        "purC" => "ade_c", # "adenine (IMP precursor)",
        "purD" => "ade_c", # "adenine (IMP precursor)",
        "purF" => "ade_c", # "adenine (IMP precursor)",
        "purL" => "ade_c", # "adenine (IMP precursor)",
        "purM" => "ade_c", # "adenine (IMP precursor)",
        "pyrB" => "orot_c", # "pyrimidine (likely orotate or uracil)",
        "pyrC" => "orot_c", # "pyrimidine (dihydroorotate? orotate?)",
        "pyrD" => "orot_c", # "pyrimidine (dihydroorotate dehydrogenase → orotate)",
        "pyrE" => "uri_c", # "uridine (orotate to UMP)",
        "pyrF" => "uri_c", # "uridine (orotidine to UMP)",
        "lysA" => "lys__L_c", # "lysine",
        "tyrA" => "tyr__L_c", # "tyrosine",
        "trpA" => "tyr__L_c", # "tryptophan",
        "trpB" => "tyr__L_c", # "tryptophan",
        "trpC" => "tyr__L_c", # "tryptophan",
        "trpD" => "tyr__L_c", # "tryptophan",
        "trpE" => "tyr__L_c", # "tryptophan",
        "ilvC" => "", # "isoleucine, leucine, or valine (one of branched-chain AAs)",
        "ilvE" => "", # "isoleucine, leucine, valine",
        "leuB" => "leu__L_c", # "leucine",
        "leuD" => "leu__L_c", # "leucine",
        "leuC" => "leu__L_c", # "leucine",
    )
end

## . -. - -.-- . . - .- .--. - < / -. .- - 
# MARK: Biomass yield
let
    global net0 = pull_net("iML1515")

    biom = extras(net0, "BIOM")
    aa_ids = [
        # Aminoacid
        "ala__L_c", "arg__L_c", "asn__L_c", "asp__L_c", "cys__L_c", "gln__L_c", "glu__L_c", "gly_c", "his__L_c", "ile__L_c", "leu__L_c", "lys__L_c", "met__L_c", "phe__L_c", "pro__L_c", "ser__L_c", "thr__L_c", "trp__L_c", "tyr__L_c", "val__L_c"
    ]


    global AA_BIOM_Y = Dict()
    for id in aa_ids
        AA_BIOM_Y[id] = stoi(net0, id, biom)
    end
    AA_BIOM_Y
end


## . -. - -.-- . . - .- .--. - < / -. .- - 
# MARK: trade yield vs biomass yield
let
    AA_BIOM_Y
    AUX_MAP
    

    aa_ids0 = String[]

    # find common ids
    for (muti, mutid) in enumerate(names(ΔN_mCherry))
        aa0 = get(AUX_MAP, mutid, "")
        isempty(aa0) && continue
        biom_y = get(AA_BIOM_Y, aa0, nothing)
        isnothing(biom_y) && continue
        @show aa0 biom_y


        push!(aa_ids0, aa0)
        
    end
    unique!(aa_ids0)

    # collect biom yield vectors
    AA_BIOM_Y_DICT = Dict()
    for aa0 in aa_ids0
        y0 = AA_BIOM_Y[aa0]
        for aa1 in aa_ids0
            aa0 == aa1 && continue
            y1 = AA_BIOM_Y[aa1]
            AA_BIOM_Y_DICT[(aa0, aa1)] = y0 / y1
        end
    end
    AA_BIOM_Y_DICT

    # collect trade yield vectos
    AA_TRADE_Y_DICT = Dict()
    for (muti0, mutid0) in enumerate(names(ΔN_mCherry))
        aa0 = get(AUX_MAP, mutid0, "")

        for (muti1, mutid1) in enumerate(names(ΔN_mCherry))
            mutid1 == mutid0 && continue
            aa1 = get(AUX_MAP, mutid1, "")
            
            n0 = ΔN_mCherry[muti0, muti1]
            n1 = ΔN_Venus[muti0, muti1]
            AA_TRADE_Y_DICT[(aa0, aa1)] = n0 / n1
        end
        # n0 = ΔN_mCherry[]
        # n1 = ΔN_Venus[]
    end
    AA_TRADE_Y_DICT


    f = Figure()
    ax = Axis(f[1,1];
        xlabel = "ΔX_i / ΔX_k (experiments)",
        ylabel = "Y_i / Y_k (biomass formulation)",
    )

    for aa0 in aa_ids0
        ys = Float64[]
        for aa1 in aa_ids0
            aa0 == aa1 && continue
            trade_y = get(AA_TRADE_Y_DICT, (aa0, aa1), nothing)
            isfinite(trade_y) || continue
            biom_y = get(AA_BIOM_Y_DICT, (aa0, aa1), nothing)
            isfinite(biom_y) || continue
            @show trade_y, biom_y
            scatter!(ax, log(trade_y), log(biom_y))
        end
    end
    f

end




## . -. - -.-- . . - .- .--. - < / -. .- - 
# MARK: gram balance
let 

    # A = sqrt.((1 .- F1_mCherry ./ F1_Venus).^2)
    # A = F1_mCherry ./ (F1_Venus .+ F1_mCherry)
    A = log10.(ΔN_mCherry ./ ΔN_Venus)
    # A = ΔN_mCherry ./ (ΔN_mCherry .+ ΔN_Venus)
    A = Matrix(A)
    A[diagind(A)] .= NaN
    # return A

    data = collect.(eachrow(A))

    # Plot
    fig = Figure(;size=(600,500))
    ax = Axis(fig[1, 1];
        limits = (nothing, nothing, -4, 5),
        xlabel = "strain k", 
        # ylabel = "growth fraction"
        ylabel = "log ΔX_i / ΔX_k"
    )

    # This works: data is Vector of Vectors
    for (i, v) in enumerate(data)
        # violin!(ax, [i], v)
        lines!(ax, eachindex(v), sort(v); 
            alpha = 0.85, 
            linewidth = 5
        )
        scatter!(ax, eachindex(v), sort(v); 
            alpha = 1.0, 
        )
        # lines!(ax, eachindex(v), v; alpha = 0.8)
    end

    fig

end



## . -. - -.-- . . - .- .--. - < / -. .- - 
function _trig_vec(A)
    M, N = size(A)
    v = []
    for i in 1:M, j in 1:N
        i < j || continue
        push!(v, A[i,j])
    end
    return v
end

## . -. - -.-- . . - .- .--. - < / -. .- - 
# Reproducing data
let
    S = 0.0:0.7:1e3
    iter = zip(
        _trig_vec(F1_mCherry), 
        _trig_vec(F1_Venus)
    )
    h0 = Dict()
    for (mc, ve) in iter
        dmc = descretize(S, mc)
        get!(h0, dmc, Float64[])
        push!(h0[dmc], ve)
    end

    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (x, dat) in h0
        push!(xs, x)
        push!(ys, mean(dat))
        push!(yerrs, std(dat))
    end

    f = Figure()
    ax = Axis(f[1,1])
    scatter!(ax, xs, ys)
    # errorbars!(ax, xs, ys, yerrs)
    f
end






## . -. - -.-- . . - .- .--. - < / -. .- - 
_vec(df) = vec(Matrix(df))
function _utrig(A)
    m, n = size(A)
    result = Vector{eltype(A)}()
    for i in 1:m
        for j in i:n  # upper triangle includes j ≥ i
            push!(result, A[i, j])
        end
    end
    return result
end


## . -. - -.-- . . - .- .--. - < / -. .- - 
let
    f = Figure()
    ax = Axis(f[1,1])

    
    scatter!(_utrig(ΔN_mCherry), _utrig(ΔN_Venus))
    f
end


## . -. - -.-- . . - .- .--. - < / -. .- - 
let
    datasheet = "mCherry_Counts"
    df1 = load_sheet(datasheet)

    datasheet = "Venus_Counts"
    df2 = load_sheet(datasheet)
    
end

## . -. - -.-- . . - .- .--. - < / -. .- - 
# MARK: Idex vs value
let
    datasheets = [
        "Day3_Rep2", "Day2_Rep1", "mCherry_Counts", "Day1_Rep1", "Day1_Rep2", "Venus_Counts", "ReadMe", "Day3_Rep1", "Day4_Rep1", "Day2_Rep2", "Day4_Rep2"
    ]

    datasheet1 = "mCherry_Counts"
    df1 = load_sheet(datasheet1)
    datasheet1 = "Venus_Counts"
    df3 = load_sheet(datasheet1)

    datasheet2 = "Day3_Rep2"
    df2 = load_sheet(datasheet2)

    xs0 = vec(Matrix(df1)) .+ vec(Matrix(df3))
    ys0 = vec(Matrix(df2))
    
    f = Figure()
    ax = Axis(f[1,1]; 
        # limits = (0, 1500, nothing, nothing),
        xlabel = "index",
        ylabel = "value",
    )

    sidx = sortperm(xs0)
    sidx = sidx
    xs1 = xs0[sidx]
    ys1 = ys0[sidx]
    scatter!(ax, ys1 ./ maximum(ys0); 
        label = datasheet1
    )
    scatter!(ax, xs1 ./ maximum(xs0); 
        label = datasheet2
    )
    f

end

## . -. - -.-- . . - .- .--. - < / -. .- - 
# MARK: minimal leakage
let
    datasheets = [
        "Day3_Rep2", "Day2_Rep1", "mCherry_Counts", "Day1_Rep1", "Day1_Rep2", "Venus_Counts", "ReadMe", "Day3_Rep1", "Day4_Rep1", "Day2_Rep2", "Day4_Rep2"
    ]

    
    datasheet1 = "Day1_Rep1"
    df1 = load_sheet(datasheet1)
    datasheet2 = "mCherry_Counts"
    df2 = load_sheet(datasheet2)
    datasheet2 = "Venus_Counts"
    df3 = load_sheet(datasheet2)
    sidx = nothing
    
    f = Figure()
    ax = Axis(f[1,1]; 
        # limits = (0, 1500, nothing, nothing),
        xlabel = "index",
        ylabel = "value",
    )

    row0 = 1
    for rowi in 1:size(df1, 1)
        rowi == row0 && continue

        v1 = Matrix(df1)[rowi, :]
        v2 = Matrix(df2)[rowi, :]
        v3 = Matrix(df3)[rowi, :]

        xs = (v2 .- v1)
        lines!(ax, xs; 
            label = datasheet1
        )
    end
    f

end

## . -. - -.-- . . - .- .--. - < / -. .- - 
let

    # Desity matrix
    datasheet = "mCherry_Counts"
    df1 = load_sheet(datasheet)
    Xi = Matrix(df1)

    datasheet = "Venus_Counts"
    df2 = load_sheet(datasheet)
    Xj = Matrix(df2)

    n = size(Xi, 1)

    model = Model(Ipopt.Optimizer)

    @variable(model, α[1:n] >= 1e-3)

    # Normalization constraint to avoid trivial zero solution
    @NLconstraint(model, sum(α[i]^2 for i in 1:n) == 1)

    # Define the objective
    @NLobjective(model, Min,
        sum(
            α[i]*Xj[i,j] - α[j]*Xi[i,j]
            for i in 1:n, j in 1:n
            if i < j
        )
    )

    optimize!(model)

    # Retrieve results
    α_opt = value.(α)
    obj_value = objective_value(model)

    println("Minimum residual: ", obj_value)
    println("Approximate solution α: ", α_opt)
end

## . -. - -.-- . . - .- .--. - < / -. .- - 