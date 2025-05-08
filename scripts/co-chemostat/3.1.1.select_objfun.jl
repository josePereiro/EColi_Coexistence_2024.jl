@cli objfun_id = "ΔMet.egoist.max.yield"

# MARK: ΔMet.egoist.max.yield
if objfun_id == "ΔMet.egoist.max.yield"
    set_linear_obj!(opm, 
        [   
            "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_glc__D_e_b",  
            "ΔIle:EX_glc__D_e_b",
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ],
        [
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
    )
# MARK: ΔIle.egoist.max.yield
elseif objfun_id == "ΔIle.egoist.max.yield"
    set_linear_obj!(opm, 
        [   
            "ΔMet:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "ΔIle:BIOMASS_Ec_iML1515_WT_75p37M_f",
            "Medium:EX_glc__D_e_b", 
            "ΔMet:EX_glc__D_e_b",  
            "ΔIle:EX_glc__D_e_b",
            "ΔMet:EX_ile__L_e_b",
            "ΔIle:EX_ile__L_e_b",
            "ΔMet:EX_met__L_e_b",
            "ΔIle:EX_met__L_e_b",
        ],
        [
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
    )
else
    error("Invalid objfun_id")
end