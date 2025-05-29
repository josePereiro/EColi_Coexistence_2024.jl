# MARK: using
using MetX
using MetXNetHub
using GLPK
using Ipopt

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: prepare_ecoli_core_opm
function prepare_ecoli_core_opm(args::Dict)

    net0 = pull_net("ecoli_core")
    lep0 = lepmodel(net0)
    
    if get(args, "test.fba", false)
        @time opm = FBAOpModel(lep0, GLPK.Optimizer)

        # max z1 + z2
        set_linear_obj!(opm, ["BIOMASS_Ecoli_core_w_GAM"], [1])
        @time optimize!(opm)

        println("-"^30)
        _opm_summary_ecoli_core(opm)
    end

    return (net0, lep0)
end

# -- .. - .-- .-. . .... -- -- -- .. ...
nothing