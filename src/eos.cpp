// Copyright 2018 @ Chun Shen

#include "eos.h"


EOS::EOS(const int eos_id_in) : eos_id(eos_id_in)  {
    if (eos_id == 0) {
        eos_ideal.initialize_eos();
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
    } else if (eos_id == 9) {
        eos_HQCD.initialize_eos();
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
    }
}

