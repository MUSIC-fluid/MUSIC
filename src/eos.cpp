// Copyright 2018 @ Chun Shen

#include "eos.h"
#include <iostream>

EOS::EOS(const int eos_id_in) : eos_id(eos_id_in)  {
    if (eos_id == 0) {
        ideal.initialize_eos();
        pressure_ptr    = &EOS::get_pressure_idealgas;
        temperature_ptr = &EOS::get_temperature_idealgas;
        entropy_ptr     = &EOS::get_entropy_idealgas;
        cs2_ptr         = &EOS::get_cs2_idealgas;
        dpde_ptr        = &EOS::get_dpde_idealgas;
        dpdrhob_ptr     = &EOS::get_dpdrhob_idealgas;
        muB_ptr         = &EOS::get_muB_idealgas;
        muS_ptr         = &EOS::get_muS_idealgas;
        muC_ptr         = &EOS::get_muC_idealgas;
        s2e_ptr         = &EOS::get_s2e_idealgas;
        T2e_ptr         = &EOS::get_T2e_idealgas;
        get_eps_max_ptr = &EOS::get_eps_max_idealgas;
        check_eos_ptr   = &EOS::check_eos_idealgas;
    } else if (eos_id == 1) {
        eosQ.initialize_eos();
        pressure_ptr    = &EOS::get_pressure_eosQ;
        temperature_ptr = &EOS::get_temperature_eosQ;
        entropy_ptr     = &EOS::get_entropy_eosQ;
        cs2_ptr         = &EOS::get_cs2_eosQ;
        dpde_ptr        = &EOS::get_dpde_eosQ;
        dpdrhob_ptr     = &EOS::get_dpdrhob_eosQ;
        muB_ptr         = &EOS::get_muB_eosQ;
        muS_ptr         = &EOS::get_muS_eosQ;
        muC_ptr         = &EOS::get_muC_eosQ;
        s2e_ptr         = &EOS::get_s2e_eosQ;
        T2e_ptr         = &EOS::get_T2e_eosQ;
        get_eps_max_ptr = &EOS::get_eps_max_eosQ;
        check_eos_ptr   = &EOS::check_eos_eosQ;
    } else if (eos_id >= 2 && eos_id <= 7) {
        s95p.initialize_eos(eos_id);
        pressure_ptr    = &EOS::get_pressure_s95p;
        temperature_ptr = &EOS::get_temperature_s95p;
        entropy_ptr     = &EOS::get_entropy_s95p;
        cs2_ptr         = &EOS::get_cs2_s95p;
        dpde_ptr        = &EOS::get_dpde_s95p;
        dpdrhob_ptr     = &EOS::get_dpdrhob_s95p;
        muB_ptr         = &EOS::get_muB_s95p;
        muS_ptr         = &EOS::get_muS_s95p;
        muC_ptr         = &EOS::get_muC_s95p;
        s2e_ptr         = &EOS::get_s2e_s95p;
        T2e_ptr         = &EOS::get_T2e_s95p;
        get_eps_max_ptr = &EOS::get_eps_max_s95p;
        check_eos_ptr   = &EOS::check_eos_s95p;
    } else if (eos_id == 8) {
        WB.initialize_eos();
        pressure_ptr    = &EOS::get_pressure_WB;
        temperature_ptr = &EOS::get_temperature_WB;
        entropy_ptr     = &EOS::get_entropy_WB;
        cs2_ptr         = &EOS::get_cs2_WB;
        dpde_ptr        = &EOS::get_dpde_WB;
        dpdrhob_ptr     = &EOS::get_dpdrhob_WB;
        muB_ptr         = &EOS::get_muB_WB;
        muS_ptr         = &EOS::get_muS_WB;
        muC_ptr         = &EOS::get_muC_WB;
        s2e_ptr         = &EOS::get_s2e_WB;
        T2e_ptr         = &EOS::get_T2e_WB;
        get_eps_max_ptr = &EOS::get_eps_max_WB;
        check_eos_ptr   = &EOS::check_eos_WB;
    } else if (eos_id == 9) {
        hotQCD.initialize_eos();
        pressure_ptr    = &EOS::get_pressure_hotQCD;
        temperature_ptr = &EOS::get_temperature_hotQCD;
        entropy_ptr     = &EOS::get_entropy_hotQCD;
        cs2_ptr         = &EOS::get_cs2_hotQCD;
        dpde_ptr        = &EOS::get_dpde_hotQCD;
        dpdrhob_ptr     = &EOS::get_dpdrhob_hotQCD;
        muB_ptr         = &EOS::get_muB_hotQCD;
        muS_ptr         = &EOS::get_muS_hotQCD;
        muC_ptr         = &EOS::get_muC_hotQCD;
        s2e_ptr         = &EOS::get_s2e_hotQCD;
        T2e_ptr         = &EOS::get_T2e_hotQCD;
        get_eps_max_ptr = &EOS::get_eps_max_hotQCD;
        check_eos_ptr   = &EOS::check_eos_hotQCD;
    } else if (eos_id >= 10 && eos_id <= 14) {
        neos.initialize_eos(eos_id);
        pressure_ptr    = &EOS::get_pressure_neos;
        temperature_ptr = &EOS::get_temperature_neos;
        entropy_ptr     = &EOS::get_entropy_neos;
        cs2_ptr         = &EOS::get_cs2_neos;
        dpde_ptr        = &EOS::get_dpde_neos;
        dpdrhob_ptr     = &EOS::get_dpdrhob_neos;
        muB_ptr         = &EOS::get_muB_neos;
        muS_ptr         = &EOS::get_muS_neos;
        muC_ptr         = &EOS::get_muC_neos;
        s2e_ptr         = &EOS::get_s2e_neos;
        T2e_ptr         = &EOS::get_T2e_neos;
        get_eps_max_ptr = &EOS::get_eps_max_neos;
        check_eos_ptr   = &EOS::check_eos_neos;
    } else if (eos_id == 17) {
        best.initialize_eos();
        pressure_ptr    = &EOS::get_pressure_best;
        temperature_ptr = &EOS::get_temperature_best;
        entropy_ptr     = &EOS::get_entropy_best;
        cs2_ptr         = &EOS::get_cs2_best;
        dpde_ptr        = &EOS::get_dpde_best;
        dpdrhob_ptr     = &EOS::get_dpdrhob_best;
        muB_ptr         = &EOS::get_muB_best;
        muS_ptr         = &EOS::get_muS_best;
        muC_ptr         = &EOS::get_muC_best;
        s2e_ptr         = &EOS::get_s2e_best;
        T2e_ptr         = &EOS::get_T2e_best;
        get_eps_max_ptr = &EOS::get_eps_max_best;
        check_eos_ptr   = &EOS::check_eos_best;
    } else {
        std::cout << "No EOS for eos_id = " << std::endl;
        exit(1);
    }
}

