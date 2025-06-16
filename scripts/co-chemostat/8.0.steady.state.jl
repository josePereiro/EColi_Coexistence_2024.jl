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
    using Statistics
    using Dates
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
# MARK: include
include("0.0.project.jl")
include("0.88.scoperias.utils.jl")
include("0.99.utils.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
let

    stst_bb = blobbatch!(B, "steady.state")
    
    # reset!
    rm(stst_bb)
    empty!(stst_bb)

    foreach_batch(B;
        n_tasks = 1
    ) do dyn_bb
        hasframe(dyn_bb, "dyn.sol") || return
        niters = dyn_bb[1]["context", "niters"]
        niters == 15000 || return

        for dyn_b in dyn_bb
            
            stst_b = blob!(stst_bb, dyn_b.uuid)

            # context
            stst_b["context", "."] = blobyref(dyn_b)
            
            # compute steady state
            dynsol = dyn_b["dyn.sol", "."]
            stst_b["stst.sol", "."] = ststsol = Dict()
            for (dyn_key, sol) in dynsol
                stst_key = replace(dyn_key, ".ts" => ".stst")
                
                ststsol[stst_key] = Dict(
                    "mean" => mean(last(sol, 500)),
                    "std" => std(last(sol, 500)),
                )
            end
        end

        @show niters
    end

    serialize!(stst_bb)
    
end