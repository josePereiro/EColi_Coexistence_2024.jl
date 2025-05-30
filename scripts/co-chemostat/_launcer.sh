#! /bin/bash

# RUN SIM
args=(
    --project
    "/Users/Pereiro/.julia/dev/EColi_Coexistence_2024/scripts/co-chemostat/co-chemostat-dfba-2strain-PAM_model.jl"
    prepare_opm="1.3.prepare.iML1515.2strains.opm.jl"
    objfun_id="ΔMet.egoist.max.yield"
    glc_Km=1.5                  # [mM] 
    glc_Vmax=8.9                # [mmol/g/h] Match PAM u_max
    dt_min=5e-3                 # [h]
    dt_max=5e-2                 # [h]
    X_max=10.0                  # [g/L]
    X_min=1e-3                  # [g/L]
    c_glc__D_e_t0=10.0
    s_glc__D_e_t0=10.0
    ΔIle_X_t0=0.3
    ΔMet_X_t0=0.3
    dt_t0=1e-1                  # [h]
    D_t0=0.1                   # [1/h]
    niters=800
)
julia "${args[@]}"
