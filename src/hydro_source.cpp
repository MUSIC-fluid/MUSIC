// Copyright 2019 Chun Shen

#include "data.h"
#include "hydro_source.h"
#include "hydro_source_strings.h"
#include "hydro_source_ampt.h"

#include <memory>

HydroSource::HydroSource(const InitData &DATA_in) : DATA(DATA_in) {
    if (DATA.Initial_profile == 13) {  // MC-Glauber-LEXUS
        hydro_source_ptr = std::unique_ptr<HydroSourceStrings> (
                                            new HydroSourceStrings (DATA));
    } else if (DATA.Initial_profile == 30) {  // AMPT
        hydro_source_ptr = std::unique_ptr<HydroSourceAMPT> (
                                            new HydroSourceStrings (DATA));
    }
}
