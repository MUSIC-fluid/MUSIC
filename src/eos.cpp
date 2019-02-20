// Copyright 2018 @ Chun Shen

#include "eos.h"
#include "eos_idealgas.h"
#include "eos_EOSQ.h"
#include "eos_s95p.h"
#include "eos_WB.h"
#include "eos_hotQCD.h"
#include "eos_best.h"
#include "eos_neos.h"
#include <iostream>

EOS::EOS(const int eos_id_in) : eos_id(eos_id_in)  {
    if (eos_id == 0) {
        eos_ptr = std::make_shared<EOS_idealgas> ();
        eos_ptr->initialize_eos();
    } else if (eos_id == 1) {
        eos_ptr = std::make_shared<EOS_eosQ> ();
        eos_ptr->initialize_eos();
    } else if (eos_id >= 2 && eos_id <= 7) {
        eos_ptr = std::make_shared<EOS_s95p> ();
        eos_ptr->initialize_eos(eos_id);
    } else if (eos_id == 8) {
        eos_ptr = std::make_shared<EOS_WB> ();
        eos_ptr->initialize_eos();
    } else if (eos_id == 9) {
        eos_ptr = std::make_shared<EOS_hotQCD> ();
        eos_ptr->initialize_eos();
    } else if (eos_id >= 10 && eos_id <= 14) {
        eos_ptr = std::make_shared<EOS_neos> ();
        eos_ptr->initialize_eos(eos_id);
    } else if (eos_id == 17) {
        eos_ptr = std::make_shared<EOS_BEST> ();
        eos_ptr->initialize_eos();
    } else {
        std::cout << "No EOS for eos_id = " << std::endl;
        exit(1);
    }
}

