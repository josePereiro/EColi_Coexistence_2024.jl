module EColi_Coexistence_2024

    export PaperSON_dir
    function PaperSON_dir()
        path = joinpath(
            pkgdir(EColi_Coexistence_2024), 
            "data", "PaperSON"
        )
        !isdir(path) && error("PaperSON folder missing... see https://github.com/MetabolicXploration/PaperSON")
        return path
    end

end
