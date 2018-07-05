// Copyright 2018 @ Chun Shen

#include "eos.h"


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
        s2e_ptr         = &EOS::get_s2e_idealgas;
        get_eps_max_ptr = &EOS::get_eps_max_idealgas;
        check_eos_ptr   = &EOS::check_eos_idealgas;
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
        s2e_ptr         = &EOS::get_s2e_s95p;
        get_eps_max_ptr = &EOS::get_eps_max_s95p;
        check_eos_ptr   = &EOS::check_eos_s95p;
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
        s2e_ptr         = &EOS::get_s2e_hotQCD;
        get_eps_max_ptr = &EOS::get_eps_max_hotQCD;
        check_eos_ptr   = &EOS::check_eos_hotQCD;
    } else if (eos_id == 10) {
        neos2.initialize_eos();
        pressure_ptr    = &EOS::get_pressure_neos2;
        temperature_ptr = &EOS::get_temperature_neos2;
        entropy_ptr     = &EOS::get_entropy_neos2;
        cs2_ptr         = &EOS::get_cs2_neos2;
        dpde_ptr        = &EOS::get_dpde_neos2;
        dpdrhob_ptr     = &EOS::get_dpdrhob_neos2;
        muB_ptr         = &EOS::get_muB_neos2;
        muS_ptr         = &EOS::get_muS_neos2;
        s2e_ptr         = &EOS::get_s2e_neos2;
        get_eps_max_ptr = &EOS::get_eps_max_neos2;
        check_eos_ptr   = &EOS::check_eos_neos2;
    } else {
        music_message << "No EOS for eos_id = " << eos_id;
        music_message.flush("error");
        exit(1);
    }
}

