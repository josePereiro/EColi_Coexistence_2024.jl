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

#=
    Register all glucose substitutes
=#

# -- .. - .-- .-. . .... -- -- -- .. ...
include("0.0.project.jl")
include("0.99.utils.jl")


## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Plot 1
let
    sources = collect(keys(carbon_sources))
    equivs = [carbon_sources[source]["ten.mmol.glc.max.equiv"] for source in sources]
    sidx = sortperm(equivs; rev = true)
    equivs = equivs[sidx]
    sources = sources[sidx]

    f = Figure()
    ax = Axis(f[1,1])
    scatter!(ax, eachindex(equivs), equivs)
    f
end