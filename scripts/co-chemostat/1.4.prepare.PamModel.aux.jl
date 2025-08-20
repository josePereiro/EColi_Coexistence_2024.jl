# MARK: using
using MetX
using MetXNetHub
using GLPK
using Ipopt
using JSON
using SparseArrays


## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: PAM Model
# u_ider = "EX_glc__D_e_b", 
# u_max = 8.9
function prepare_PamModel_alter_single_aux_opm(args::Dict)

    net1 = pull_net(
        "PAMModel-alterProteomeRegulationPatterns2021", 
        Dict(
            "u_ider" => args["u_ider"],
            "u_max" => args["u_max"],
        )
    )
    # bounds!(net0, "METabcpp_b", 0.0, 1000.0)
    # bounds!(net0, "METabcpp_f", 0.0, 1000.0)
    
    aux = args["aux"]
    if aux == "ΔMet"
        # HSST_f
        # (-1.0) hom__L_c + (-1.0) succoa_c ==> (1.0) coa_c + (1.0) suchms_c
        extras!(net1, "ΔMet_KO", "HSST_f")
        # close KOs
        bounds!(net1, "HSST_f", 0.0, 0.0)
        bounds!(net1, "HSST_b", 0.0, 0.0)
    elseif aux == "ΔIle"
        # THRD_L_f
        # (-1.0) thr__L_c ==> (1.0) 2obut_c + (1.0) nh4_c
        extras!(net1, "ΔIle_KO", "ΔIle:THRD_L_f")
        # close KOs
        bounds!(net1, "THRD_L_f", 0.0, 0.0)
        bounds!(net1, "THRD_L_b", 0.0, 0.0)
    else
        error("invalid aux '$(aux)'")
    end

    # But, close carbon sources
    # see ...carbon-equiv.jl
    jsonpath = _procdir(["carbon-equiv"], "iML1515.carbon_sources.json")
    carbon_sources = JSON.parsefile(jsonpath)
    for ex0 in keys(carbon_sources)
        ex1 = "$(ex0)_b" # backward are the imports
        ub!(net1, "$ex1", 0.0)
    end

    # block reverse biomass
    bounds!(net1, "BIOMASS_Ec_iML1515_WT_75p37M_b", 0.0, 0.0)
    
    # Medium
    # --------------------
    # MARK: EX_glc__D_e_f
    # (-1.0) glc__D_e ==>
    # EX_glc__D_e_b
    # ==> (1.0) glc__D_e
    bounds!(net1, "EX_glc__D_e_f", 0.0, 0.0)
    bounds!(net1, "EX_glc__D_e_b", 0.0, 1000.0)

    # --------------------
    # MARK: EX_met__L_e_f
    # (-1.0) met__L_e ==>
    # EX_met__L_e_b
    # ==> (1.0) met__L_e
    bounds!(net1, "EX_met__L_e_f", 0.0, 1000.0)
    bounds!(net1, "EX_met__L_e_b", 0.0, 0.0)

    # --------------------
    # MARK: EX_ile__L_e_f
    # (-1.0) ile__L_e ==> 
    # EX_ile__L_e_b
    # ==> (1.0) ile__L_e
    bounds!(net1, "EX_ile__L_e_f", 0.0, 1000.0)
    bounds!(net1, "EX_ile__L_e_b", 0.0, 0.0)

    # -- .. - .-- .-. . .... -- -- -- .. ...
    # MARK: lepmodel
    lep0 = lepmodel(net1)

    # -- .. - .-- .-. . .... -- -- -- .. ...
    # MARK: Test fba
    if get(args, "test.fba", false)
        @time opm = FBAOpModel(lep0, Clp.Optimizer)

        if aux == "ΔIle"
            bounds!(opm, "EX_ile__L_e_b", 0.0, 1000.0)
        elseif aux == "ΔMet"
            bounds!(opm, "EX_met__L_e_b", 0.0, 1000.0)
        end

        # max z1 + z2
        set_linear_obj!(opm, 
            ["BIOMASS_Ec_iML1515_WT_75p37M_f"], 
            [1]
        )
        @time optimize!(opm)

        println("-"^30)
        _opm_summary_PamModel_alter_single_aux(opm)
    end

    return (net1, lep0)
end