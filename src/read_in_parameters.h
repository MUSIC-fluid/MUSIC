#ifndef SRC_READ_IN_PARAMETERS_H_
#define SRC_READ_IN_PARAMETERS_H_

#include <string>

#include "data.h"
#include "util.h"
#include "emoji.h"
#include "pretty_ostream.h"

//! This class handles read in parameters
namespace ReadInParameters {
    InitData read_in_parameters(std::string input_file);
    void check_parameters(InitData &parameter_list, std::string input_file);
}

#endif  // SRC_READ_IN_PARAMETERS_H_
