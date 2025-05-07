# MARK: using
@time begin
    using MetX
    using MetX.MetXNetHub
    using MetX.GLPK
    using MetX.Ipopt
    using JSON
    using SparseArrays
end

# -- .. - .-- .-. . .... -- -- -- .. ...
include("0.utils.jl")

# -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: PAM Model
u_ider = "EX_glc__D_e_b"
u_max = 8.9
net0 = pull_net(
    "PAMModel-alterProteomeRegulationPatterns2021", 
    Dict(
        "u_ider" => u_ider,
        "u_max" => u_max,
    )
)
# bounds!(net0, "METabcpp_b", 0.0, 1000.0)
# bounds!(net0, "METabcpp_f", 0.0, 1000.0)

net1 = double_model(net0, "ΔMet", "ΔIle", "Medium")

# HSST_f
# (-1.0) hom__L_c + (-1.0) succoa_c ==> (1.0) coa_c + (1.0) suchms_c
extras!(net1, "ΔMet_KO", "ΔMet:HSST_f")
# THRD_L_f
# (-1.0) thr__L_c ==> (1.0) 2obut_c + (1.0) nh4_c
extras!(net1, "ΔIle_KO", "ΔIle:THRD_L_f")

# close KOs
bounds!(net1, "ΔMet:HSST_f", 0.0, 0.0)
bounds!(net1, "ΔMet:HSST_b", 0.0, 0.0)
bounds!(net1, "ΔIle:THRD_L_f", 0.0, 0.0)
bounds!(net1, "ΔIle:THRD_L_b", 0.0, 0.0)

# Open inter-cells exchanges
exchs0 = extras(net1, "EX_RXNS0")
for ex0 in exchs0
    bounds!(net1, "ΔMet:$ex0", 0.0, 1000.0)
    bounds!(net1, "ΔIle:$ex0", 0.0, 1000.0)
end

# But, close carbon sources
jsonpath = joinpath(@__DIR__, "carbon_sources.json")
carbon_sources = JSON.parsefile(jsonpath)
for ex0 in keys(carbon_sources)
    ex1 = "$(ex0)_b" # backward are the imports
    ub!(net1, "ΔMet:$ex1", 0.0)
    ub!(net1, "ΔIle:$ex1", 0.0)
end

# block reverse biomass
bounds!(net1, "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_b", 0.0, 0.0)
bounds!(net1, "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_b", 0.0, 0.0)

# Medium
# --------------------
# MARK: X:EX_glc__D_e_f
# (-1.0) X:glc__D_e ==> (1.0) Medium:glc__D_e
# X:EX_glc__D_e_b
# (-1.0) Medium:glc__D_e ==> (1.0) X:glc__D_e
bounds!(net1, "ΔMet:EX_glc__D_e_f", 0.0, 0.0)
bounds!(net1, "ΔMet:EX_glc__D_e_b", 0.0, 1000.0)
bounds!(net1, "ΔIle:EX_glc__D_e_f", 0.0, 0.0)
bounds!(net1, "ΔIle:EX_glc__D_e_b", 0.0, 1000.0)

# Medium:EX_glc__D_e_f
# (-1.0) Medium:glc__D_e ==>
# Medium:EX_glc__D_e_b
# ==> (1.0) Medium:glc__D_e
bounds!(net1, "Medium:EX_glc__D_e_f", 0.0, 0.0)
bounds!(net1, "Medium:EX_glc__D_e_b", 0.0, 1000.0)

# --------------------
# MARK: X:EX_met__L_e_f
# (-1.0) X:met__L_e ==> (1.0) Medium:met__L_e
# X:EX_met__L_e_b
# (-1.0) Medium:met__L_e ==> (1.0) X:met__L_e
bounds!(net1, "ΔMet:EX_met__L_e_f", 0.0, 1000.0)
bounds!(net1, "ΔMet:EX_met__L_e_b", 0.0, 1000.0)
bounds!(net1, "ΔIle:EX_met__L_e_f", 0.0, 1000.0)
bounds!(net1, "ΔIle:EX_met__L_e_b", 0.0, 1000.0)

# Medium:EX_met__L_e_f
#  (-1.0) Medium:met__L_e ==> 
# Medium:EX_met__L_e_b 
# ==> (1.0) Medium:met__L_e
bounds!(net1, "Medium:EX_met__L_e_f", 0.0, 1000.0)
bounds!(net1, "Medium:EX_met__L_e_b", 0.0, 0.0)

# --------------------
# MARK: X:EX_ile__L_e_f
# (-1.0) X:ile__L_e ==> (1.0) Medium:ile__L_e
# X:EX_ile__L_e_b
# (-1.0) Medium:ile__L_e ==> (1.0) X:ile__L_e
bounds!(net1, "ΔMet:EX_ile__L_e_f", 0.0, 1000.0)
bounds!(net1, "ΔMet:EX_ile__L_e_b", 0.0, 1000.0)
bounds!(net1, "ΔIle:EX_ile__L_e_f", 0.0, 1000.0)
bounds!(net1, "ΔIle:EX_ile__L_e_b", 0.0, 1000.0)

# Medium:EX_ile__L_e_f
# (-1.0) Medium:ile__L_e ==>
# Medium:EX_ile__L_e_b
# ==> (1.0) Medium:ile__L_e
bounds!(net1, "Medium:EX_ile__L_e_f", 0.0, 1000.0)
bounds!(net1, "Medium:EX_ile__L_e_b", 0.0, 0.0)

# -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: Test fba
@time opm = FBAOpModel(lepmodel(net1), Clp.Optimizer)
# @time opm = FBAOpModel(lepmodel(net1), Ipopt.Optimizer)

# max z1 + z2
set_linear_obj!(opm, 
    [
        "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f", 
        "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
    ], 
    [1, 1]
)
@time optimize!(opm)

println("-"^30)
_opm_summary_PamModel_alter(opm)
