// Copyright 2018 @ Chun Shen

#include "eos.h"

#include <iostream>
#include <memory>

#include "eos_1DGenerator.h"
#include "eos_EOSQ.h"
#include "eos_UH.h"
#include "eos_WB.h"
#include "eos_best.h"
#include "eos_hotQCD.h"
#include "eos_idealgas.h"
#include "eos_neos.h"
#include "eos_s95p.h"

// NOTE: EOS_chem (eos_id 16, see eos_chem.h) is intentionally NOT wired up
// here. Y_q enters the EOS interface through the third argument added in
// eos_base.h, which defaults to 1 (full chemical equilibrium). The generic
// call sites in the evolution code (reconst.cpp, dissipative.cpp,
// eos_base.cpp, and others) omit that argument, so if EOS_chem were
// selectable as the primary EOS_to_use they would always evaluate its
// Y_q = 1 slice -- which is identical to the lattice EOS they already use,
// only slower. The physically relevant Y_q of a fluid cell is not a
// conserved density that flows through those call sites; it is derived
// from the evolved pi_b_chem field (see the Y_q calculation in
// grid_info.cpp) and must be passed explicitly. Construct EOS_chem
// directly wherever a real Y_q value is available instead.
EOS::EOS(const int eos_id_in) : eos_id(eos_id_in) {
    if (eos_id == 0) {
        eos_ptr = std::unique_ptr<EOS_idealgas>(new EOS_idealgas());
    } else if (eos_id == 1) {
        eos_ptr = std::unique_ptr<EOS_eosQ>(new EOS_eosQ());
    } else if (eos_id >= 2 && eos_id <= 7) {
        eos_ptr = std::unique_ptr<EOS_s95p>(new EOS_s95p(eos_id));
    } else if (eos_id == 8) {
        eos_ptr = std::unique_ptr<EOS_WB>(new EOS_WB());
    } else if (eos_id == 9 || eos_id == 91) {
        eos_ptr = std::unique_ptr<EOS_hotQCD>(new EOS_hotQCD(eos_id));
    } else if (eos_id >= 10 && eos_id <= 15) {
        eos_ptr = std::unique_ptr<EOS_neos>(new EOS_neos(eos_id));
    } else if (eos_id == 17) {
        eos_ptr = std::unique_ptr<EOS_BEST>(new EOS_BEST());
    } else if (eos_id == 19) {
        eos_ptr = std::unique_ptr<EOS_UH>(new EOS_UH());
    } else if (eos_id == 42) {
        eos_ptr = std::unique_ptr<EOS_1DGenerator>(new EOS_1DGenerator(eos_id));
    } else if (eos_id == 16) {
        std::cout << "EOS_to_use 16 (EOS_chem) cannot be used as the "
                   << "primary EOS -- see the NOTE above EOS::EOS in "
                   << "eos.cpp for why. Construct EOS_chem directly instead."
                   << std::endl;
        exit(1);
    } else {
        std::cout << "No EOS for eos_id = " << std::endl;
        exit(1);
    }
    eos_ptr->initialize_eos();
}
