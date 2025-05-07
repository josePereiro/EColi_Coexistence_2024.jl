# MARK: using
@time begin
    using MetX
    using MetX.MetXNetHub
    using MetX.GLPK
    using MetX.Ipopt
    using JSON
end

# -- .. - .-- .-. . .... -- -- -- .. ...
include("0.utils.jl")

## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Model
net0 = pull_net("ecoli_core")

@time opm = FBAOpModel(lepmodel(net0), GLPK.Optimizer)

# max z1 + z2
set_linear_obj!(opm, ["BIOMASS_Ecoli_core_w_GAM"], [1])
@time optimize!(opm)

println("-"^30)
_opm_summary_ecoli_core(opm)

# -- .. - .-- .-. . .... -- -- -- .. ...
nothing