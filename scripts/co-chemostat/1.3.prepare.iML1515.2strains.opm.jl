# MARK: using
using MetX
using MetXNetHub
using GLPK
using Ipopt
using JSON
using SparseArrays

# -- .. - .-- .-. . .... -- -- -- .. ...
function prepare_iML1515_2strains_opm(args::Dict)

    # MARK: PAM Model
    net0 = pull_net("iML1515")
    # bounds!(net0, "METabcpp_b", 0.0, 1000.0)
    # bounds!(net0, "METabcpp_f", 0.0, 1000.0)

    u_max = 8.9 # nutrient max intake rate

    # Extract
    S0, b0, lb0, ub0, mets0, rxns0 = net0.S, net0.b, net0.lb, net0.ub, net0.mets, net0.rxns

    # Split
    rxns1 = [
        [string(rxn, "_f") for rxn in rxns0]; 
        [string(rxn, "_b") for rxn in rxns0]; 
    ]
    mets1 = mets0
    S1 = [
        S0  -S0
    ]
    M, N = size(S1)
    b1 = b0
    c1 = zeros(N)
    lb1 = zeros(N)
    ub1 = [ub0; abs.(lb0)]

    # MARK: MetNet
    net1 = MetNet(;S=S1, b=b1, lb=lb1, ub=ub1, c=c1, rxns=rxns1, mets=mets1)

    extras!(net1, "BIOM", "BIOMASS_Ec_iML1515_WT_75p37M_f")
    extras!(net1, "EX_GLC", "EX_glc__D_e_b")
    extras!(net1, "EX_NH4", "EX_nh4_e_b")
    extras!(net1, "EX_GLU", "EX_glu__L_e_b")
    extras!(net1, "EX_O2", "EX_o2_e_b")
    extras!(net1, "EX_CO2", "EX_co2_e_b")
    extras!(net1, "ATPM", "ATPM")

    net1 = double_model(net1, "ΔMet", "ΔIle", "Medium")


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
    # see ...carbon-equiv.jl
    jsonpath = _procdir(["carbon-equiv"], "iML1515.carbon_sources.json")
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
    bounds!(net1, "ΔMet:EX_glc__D_e_b", 0.0, u_max)
    bounds!(net1, "ΔIle:EX_glc__D_e_f", 0.0, 0.0)
    bounds!(net1, "ΔIle:EX_glc__D_e_b", 0.0, u_max)

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
    # MARK: lep
    lep0 = lepmodel(net1)

    # -- .. - .-- .-. . .... -- -- -- .. ...
    # MARK: Test fba
    if get(args, "test.fba", false)
        @time opm = FBAOpModel(lep0, Clp.Optimizer)
        # @time opm = FBAOpModel(lep0, Ipopt.Optimizer)

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
    end

    return (net1, lep0)

end