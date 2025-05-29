# MARK: using
@time begin
    using CSV
    using MetX
    using JSON
    using Clp
    using GLPK
    using Ipopt
    using MetXNetHub
    using DataFrames
    using SparseArrays
    using Serialization
    using EColi_Coexistence_2024
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("0.0.project.jl")
include("0.99.utils.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("1.0.prepare.lepmodel.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
let
    # @st_flush!
    lepmodel_id = "core_ecoli"
    net0, opm = prepare_lepmodel(lepmodel_id)
    nothing
end