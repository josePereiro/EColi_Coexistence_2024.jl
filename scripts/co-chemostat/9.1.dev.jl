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
# About
# - 


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("0.0.project.jl")
include("0.99.utils.jl")
include("9.99.utils.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: make files
let
    bb =  B[r"glc.aa.phy"]

    # filter control
    for (aux_id, nonlim_arg) in Iterators.product(
        ["ΔMet", "ΔIle"],
        ["glc_U", "aa_U"]
    )

        # U_lines
        U_LINES = U_lines(aux_id)

        # Env
        ENV = U_envelop(U_LINES, nonlim_arg)

        # write envelop
        fn = _procdir("growth-curve", (;aux_id, nonlim_arg), ".csv")
        @show basename(fn)
        open(fn, "w") do io
            println(io, nonlim_arg, ", " , "biom")
            xs = sort(collect(keys(ENV)))
            for x in xs
                b = ENV[x]
                println(io, x, ",", b)
            end
        end
    end

end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
let
    # U_lines
    aux_id = "ΔIle"
    nonlim_arg = "glc_U"
    
    @time global U_LINES = U_lines(aux_id)
end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: growth vs glc
let

    # U_lines
    aux_id = "ΔIle"
    nonlim_arg = "aa_U"
    
    @time global U_LINES = U_lines(aux_id)
    @time global ENVLOP = U_envelop(U_LINES, nonlim_arg)

    qs = ENVLOP[:xs]
    zs = ENVLOP[:ys]
    qm, qM = extrema(qs)
    @show qm, qM

    Q = qM 
    # K = 0.39 * 1e-3 # mM glc
    K = 0.39 * 1e-3 # mM glc
    @assert Q >= qM
    
    # add more
    qes = Float64[]
    zes = Float64[]
    for it in 1:5000
        xt = qm + (qM - qm) * rand()
        ye = extrapolate(xt, qs, zs)
        push!(qes, xt)
        push!(zes, ye)
    end

    # 
    ss = map(qs) do q
        s = q * K / (Q - q)
    end
    ses = map(qes) do q
        s = q * K / (Q - q)
    end

    # write envelop
    fn = _procdir("fba", (;aux_id, nonlim_arg), ".csv")
    @show basename(fn)
    open(fn, "w") do io
        println(io, 
            join(["growth", "q", "s", "yield"], ", ")
        )
        xs = sort(collect(keys(ENV)))
        for i in eachindex(zes)
            dat = [zes[i], qes[i], ses[i], zes[i] ./ qes[i]]
            println(io, join(dat, ", "))
        end
    end


    # DATA from shehataEffectNutrientConcentration1971 Fig1
    # growth vs glc
    exp_xs = [0.0815647934185284, 0.4890220252875088, 0.9916738399775272, 4.963171709158358, 9.90697539429683, 19.916077542337593]

    exp_ys = [0.16168445769943132, 0.43226748160376127, 0.5750919065344462, 0.7176572059376501, 0.7521754040910985, 0.7569141195123951]


    f = Figure()
    ax = Axis(f[1,1];
        title = aux_id,
        limits = (-1 * 1e-3, 30 * 1e-3, -0.1, 0.8),
        xlabel = "s_glc [mM]",
        ylabel = "growth [1/h]"
    )
    # scatter!(ax, ss, zs;
    #     label = "fba"
    # )
    scatter!(ax, ses, zes;
        label = "fba"
    )
    # scatter!(ax, exp_xs .* 1e-3, exp_ys;
    #     marker = :star5,
    #     markersize = 20,
    #     label = "experimental"
    # )
    # scatter!(ax, log10.(ses), zes)
    axislegend(ax; position = :rb)
    f

end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: yield(s)
let

    # U_lines
    aux_id = "ΔIle"
    nonlim_arg = "aa_U"
    
    @time global U_LINES = U_lines(aux_id)
    @time global ENVLOP = U_envelop(U_LINES, nonlim_arg)

    qs = ENVLOP[:xs]
    zs = ENVLOP[:ys]
    qm, qM = extrema(qs)
    @show qm, qM

    fn = "/Users/Pereiro/.julia/dev/EColi_Coexistence_2024/data/schulzeRelationshipSubstrateConcentration1964-table1.csv"
    exp_df = CSV.read(fn, DataFrame)    

    exp_xs = exp_df[1:end-2, "Substrate conc. (mg/l)"] ./ 180
    exp_ys = exp_df[1:end-2, "Yield factor Y"] .* 180 .* 1e-3
    @show exp_ys

    # "D (hr^-1)"
    # "Substrate conc. (mg/l)"
    # "Cell conc. (mg/l)"
    # "QO2 (mg O2/g cell wt/hr)"
    # "Yield factor Y"
    # "Glucose oxidized (%)"
    # "Rate of substrate uptake (mg/g cell wt/hr)"

    Q = qM 
    # K = 0.39 * 1e-3 # mM glc
    # K = 0.39 * 1e-3 # mM glc
    K = 0.39 # mM glc
    @assert Q >= qM
    
    # add more
    qes = Float64[]
    zes = Float64[]
    for it in 1:5000
        xt = qm + (qM - qm) * rand()
        ye = extrapolate(xt, qs, zs)
        push!(qes, xt)
        push!(zes, ye)
    end



    # 
    ss = map(qs) do q
        s = q * K / (Q - q)
    end
    ses = map(qes) do q
        s = q * K / (Q - q)
    end

    f = Figure()
    ax = Axis(f[1,1];
        title = aux_id,
        # limits = (-1 * 1e-3, 30 * 1e-3, -0.1, 0.8),
        limits = (-0.1, 1 , nothing, nothing),
        xlabel = "s_glc [mM]",
        ylabel = "growth [1/h]"
    )
    # scatter!(ax, ss, zs;
    #     label = "fba"
    # )
    scatter!(ax, ses, zes ./ qes;
        label = "fba"
    )
    scatter!(ax, exp_xs, exp_ys;
        marker = :star5,
        markersize = 20,
        label = "experimental"
    )
    # scatter!(ax, log10.(ses), zes)
    axislegend(ax; position = :rb)
    f

end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: yield(D)
let
    fn = "/Users/Pereiro/.julia/dev/EColi_Coexistence_2024/data/schulzeRelationshipSubstrateConcentration1964-table1.csv"
    exp_df = CSV.read(fn, DataFrame)    

    names(exp_df)

    # "D (hr^-1)"
    # "Substrate conc. (mg/l)"
    # "Cell conc. (mg/l)"
    # "QO2 (mg O2/g cell wt/hr)"
    # "Yield factor Y"
    # "Glucose oxidized (%)"
    # "Rate of substrate uptake (mg/g cell wt/hr)"

    # xs = exp_df[:, "Substrate conc. (mg/l)"] / 180
    exp_xs = exp_df[:, "D (hr^-1)"] 
    # [gX/gglc] / [g/mol] = []
    # 
    exp_ys = exp_df[:, "Yield factor Y"] .* 180 .* 1e-3
    # ys = exp_df[:, "D (hr^-1)"]


    # U_lines
    aux_id = "ΔIle"
    nonlim_arg = "aa_U"
    
    @time U_LINES = U_lines(aux_id)
    # flux, growth
    @time ENVLOP = U_envelop(U_LINES, nonlim_arg)
    ys, xs = ENVLOP

    f = Figure()
    ax = Axis(f[1,1];
        title = aux_id,
        xlabel = "D [1/h]",
        ylabel = "yield [g/mmol]"
    )
    scatter!(exp_xs, exp_ys;
        markersize = 20,
        label = "experiment"
    )
    scatter!(xs, xs ./ ys;  
        markersize = 20,
        label = "fba"
    )
    axislegend(ax; position = :rt)
    f
end
