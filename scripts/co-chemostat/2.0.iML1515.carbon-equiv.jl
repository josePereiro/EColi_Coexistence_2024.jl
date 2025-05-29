# MARK: using
@time begin
    using MetX
    using MetXNetHub
    using GLPK
    using Ipopt
    using CSV
    using DataFrames
    using JSON
    using Serialization
    using SparseArrays
    using Random
end

#=
    Register all glucose substitutes
=#

# -- .. - .-- .-. . .... -- -- -- .. ...
include("0.0.project.jl")
include("0.99.utils.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
let
    net0 = pull_net("iML1515")

    # Open exports
    for exch in reactions(net0)
        contains(exch, "EX_") || continue
        ub!(net0, exch, 1000.0)
    end

    ## -- .. - .-- .-. . .... -- -- -- .. ...
    # MARK: test glc
    @time opm = FBAOpModel(lepmodel(net0), GLPK.Optimizer)
    set_linear_obj!(opm, ["BIOMASS_Ec_iML1515_WT_75p37M"], [1])

    # open glc
    bounds!(opm, "EX_glc__D_e", -10.0, 1000)
    biom0 = nothing
    try
        optimize!(opm)
        biom0 = solution(opm, "BIOMASS_Ec_iML1515_WT_75p37M")
        @show biom0
    catch e
        println("Error: $e")
    end


    # fba
    bounds!(opm, "EX_glc__D_e", 0, 1000)
    try
        optimize!(opm)
        biom = solution(opm, "BIOMASS_Ec_iML1515_WT_75p37M")
        @show biom
        @assert biom < 1e-3
    catch e
        println("Error: $e")
    end


    ## -- .. - .-- .-. . .... -- -- -- .. ...
    # MARK: find equivs
    carbon_sources = Dict()

    for exch in reactions(net0)
        contains(exch, "EX_") || continue
        lb(opm, exch) == 0 || continue # already open (ignore)
        # fba
        try
            lb!(opm, exch, -10.0)
            optimize!(opm)
            biom = solution(opm, "BIOMASS_Ec_iML1515_WT_75p37M")
            @show exch
            @show biom
            biom > 0.0 || continue
            carbon_sources[exch] = Dict(
                "biom" => biom,
                "ten.mmol.glc.max.equiv" => biom / biom0,
            )
        catch e
            (e isa InterruptException) && break
            println("Error: $e")
        finally
            # close the exchange
            lb!(opm, exch, 0.0)
        end
    end

    ## -- .. - .-- .-. . .... -- -- -- .. ...
    # Save carbon_sources into a json
    datdir = _procdir("carbon-equiv")
    mkpath(datdir)
    jsonpath = joinpath(datdir, "iML1515.carbon_sources.json")
    open(jsonpath, "w") do io
        JSON.print(io, carbon_sources, 2)
    end
    
end
