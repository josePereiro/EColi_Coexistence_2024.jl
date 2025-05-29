
## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
function select_objfun(lepmodel_id, objfun_id)

    # @cli objfun_id = "ΔMet.egoist.max.yield"

    # MARK: iML1515 : ΔMet.egoist.max.yield
    if lepmodel_id in [
            "PamModel-alter.2strains", 
            "iML1515.2strains"
        ] &&
        objfun_id in [
            "ΔMet.egoist.max.yield"
        ]

        objiders = [   
            "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_glc__D_e_b",  
            "ΔIle:EX_glc__D_e_b",
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ]
        objcoes = [
                1e4,  # "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f" max ~1.0
                1e2,  # "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
                0.0,  # "Medium:EX_glc__D_e_b" 
            -1e-3,  # "ΔMet:EX_glc__D_e_b" max ~10.0
            -1e-3,  # "ΔIle:EX_glc__D_e_b"
            -1e-5,  # "ΔMet:EX_ile__L_e_b"
            -1e-5,  # "ΔIle:EX_ile__L_e_b"
            -1e-5,  # "ΔMet:EX_met__L_e_b"
            -1e-5,  # "ΔIle:EX_met__L_e_b"
        ]
        return (objiders, objcoes)

    # MARK: iML1515 : ΔIle.egoist.max.yield
    elseif lepmodel_id in [
            "PamModel-alter.2strains", 
            "iML1515.2strains"
        ] &&
        objfun_id in [
            "ΔIle.egoist.max.yield"
        ]
        objiders = [   
            "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_glc__D_e_b",  
            "ΔIle:EX_glc__D_e_b",
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ]
        objcoes = [
            1e2,  # "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f" max ~1.0
            1e4,  # "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
            0.0,  # "Medium:EX_glc__D_e_b" 
            -1e-3,  # "ΔMet:EX_glc__D_e_b" max ~10.0
            -1e-3,  # "ΔIle:EX_glc__D_e_b"
            -1e-5,  # "ΔMet:EX_ile__L_e_b"
            -1e-5,  # "ΔIle:EX_ile__L_e_b"
            -1e-5,  # "ΔMet:EX_met__L_e_b"
            -1e-5,  # "ΔIle:EX_met__L_e_b"
        ]

        return (objiders, objcoes)

    # MARK: iML1515 : ΔMet.egoist.only
    elseif lepmodel_id in [
            "PamModel-alter.2strains", 
            "iML1515.2strains"
        ] &&
        objfun_id in [
            "ΔMet.egoist.only"
        ]

        objiders = [   
            "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_glc__D_e_b",  
            "ΔIle:EX_glc__D_e_b",
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ]
        objcoes = [
                1e4,  # "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f" max ~1.0
                0.0,  # "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f"
                0.0,  # "Medium:EX_glc__D_e_b" 
            -1e-3,  # "ΔMet:EX_glc__D_e_b" max ~10.0
            -1e-3,  # "ΔIle:EX_glc__D_e_b"
            -1e-5,  # "ΔMet:EX_ile__L_e_b"
            -1e-5,  # "ΔIle:EX_ile__L_e_b"
            -1e-5,  # "ΔMet:EX_met__L_e_b"
            -1e-5,  # "ΔIle:EX_met__L_e_b"
        ]
        return (objiders, objcoes)

    # MARK: iML1515 : ΔIle.egoist.only
    elseif lepmodel_id in [
            "PamModel-alter.2strains", 
            "iML1515.2strains"
        ] &&
        objfun_id in [
            "ΔIle.egoist.only"
        ]
        
        objiders = [   
            "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f", 
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_glc__D_e_b",  
            "ΔIle:EX_glc__D_e_b",
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ]
        objcoes = [
                0.0,  # "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f" 
                1e4,  # "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f" max ~1.0
                0.0,  # "Medium:EX_glc__D_e_b" 
            -1e-3,  # "ΔMet:EX_glc__D_e_b" max ~10.0
            -1e-3,  # "ΔIle:EX_glc__D_e_b"
            -1e-5,  # "ΔMet:EX_ile__L_e_b"
            -1e-5,  # "ΔIle:EX_ile__L_e_b"
            -1e-5,  # "ΔMet:EX_met__L_e_b"
            -1e-5,  # "ΔIle:EX_met__L_e_b"
        ]
        
        return (objiders, objcoes)

    # MARK: iML1515 : non-competitive.max.yield
    elseif lepmodel_id in [
            "PamModel-alter.2strains", 
            "iML1515.2strains"
        ] &&
        objfun_id in [
            "non-competitive.max.yield"
        ]
        
        objiders = [   
            "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f", 
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_glc__D_e_b",  
            "ΔIle:EX_glc__D_e_b",
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ]
        objcoes = [
                1e4,  # "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f" 
                1e4,  # "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f" max ~1.0
                0.0,  # "Medium:EX_glc__D_e_b" 
            -1e-3,  # "ΔMet:EX_glc__D_e_b" max ~10.0
            -1e-3,  # "ΔIle:EX_glc__D_e_b"
            -1e-5,  # "ΔMet:EX_ile__L_e_b"
            -1e-5,  # "ΔIle:EX_ile__L_e_b"
            -1e-5,  # "ΔMet:EX_met__L_e_b"
            -1e-5,  # "ΔIle:EX_met__L_e_b"
        ]
        
        return (objiders, objcoes)

    # MARK: ecoli_core : max.biomass.max.yield
    elseif lepmodel_id in [
            "core_ecoli", 
        ] &&
        objfun_id in [
            "max.biomass.max.yield"
        ]

        objiders = [
            "BIOMASS_Ecoli_core_w_GAM",
            "EX_glc__D_e",
            "EX_nh4_e",
            "EX_o2_e",
            "EX_co2_e"
        ]

        objcoes = [
            1e5,
            -1e-2,
            -1e-2,
            -1e-2,
            -1e-2
        ]
        
        return (objiders, objcoes)

    else
        error("Invalid (lepmodel_id, objfun_id), ($(lepmodel_id), $(objfun_id))")
    end

end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
nothing