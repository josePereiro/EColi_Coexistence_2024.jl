## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
include("1.1.prepare.PamModel-alter.2strains.opm.jl")
include("1.2.prepare.core_ecoli.opm.jl")
include("1.3.prepare.iML1515.2strains.opm.jl")

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
function prepare_lepmodel(lepmodel_id::String, args::Dict = Dict())
    if lepmodel_id == "PamModel-alter.2strains"
        prepare_PamModel_alter_2strains_opm(args)
    elseif lepmodel_id == "core_ecoli"
        prepare_ecoli_core_opm(args)
    elseif lepmodel_id == "iML1515.2strains"
        prepare_iML1515_2strains_opm(args)
    else
        error("Unknown lepmodel_id, $(lepmodel_id)")
    end
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
nothing