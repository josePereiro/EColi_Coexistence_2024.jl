#! /bin/bash

# RUN SIM
Ds=(0.6)
dts=(1e-2, 1e-3, 1e-4)

# objfun_id0="ΔMet.egoist.max.yield"
objfun_id0="ΔIle.egoist.max.yield"
# objfun_id0="ΔIle.egoist.only"

# prepare_opm0="1.1.prepare.PamModel-alter.2strains.opm.jl"
prepare_opm0="1.3.prepare.iML1515.2strains.opm.jl"

for ((Di=0; Di<${#Ds[@]}; Di+=1)); do
for ((dti=0; dti<${#dts[@]}; dti+=1)); do
    # Get the current number and the next one
    args=(
        --project
        "/Users/Pereiro/.julia/dev/EColi_Coexistence_2024/scripts/co-chemostat/3.1.co-chemostat-dfba-2strain-PAM_model.jl"
        prepare_opm=${prepare_opm0}
        objfun_id=${objfun_id0}
        glc_Km=1.5                  # [mM] 
        glc_Vmax=8.9                # [mmol/g/h] Match PAM u_max
        dt_min=1e-4                 # [h]
        dt_max=1e-2                 # [h]
        X_max=10.0                  # [g/L]
        X_min=1e-3                  # [g/L]
        c_glc__D_e_t0=10.0
        s_glc__D_e_t0=10.0
        ΔIle_X_t0=0.3
        ΔMet_X_t0=0.3
        dt_t0=${dts[dti]}                  # [h]
        D_t0=${Ds[Di]}               # [1/h]
        niters=500
    )
    julia "${args[@]}"
done
done

