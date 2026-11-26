#ifndef SRC_READ_IN_PARAMETERS_H_
#define SRC_READ_IN_PARAMETERS_H_

#include <string>

#include "data.h"
#include "emoji.h"
#include "pretty_ostream.h"
#include "util.h"

//! This class handles read in parameters
namespace ReadInParameters {
InitData read_in_parameters(std::string input_file);
void check_parameters(InitData &parameter_list, std::string input_file);
void set_parameter(
    InitData &parameter_list, std::string parameter_name, double value);

//! the bookkeeping file of every run: music_parameters.dat
std::string parameterFileName();
//! write the parameters of the run: [input] every parameter looked up in
//! input_file with its value or default (an input file for the same run),
//! [unused] the input entries that are not read, and [run] all of DATA
void writeParameterFile(
    const InitData &DATA, const std::string &input_file,
    const std::string &filename = parameterFileName());
}  // namespace ReadInParameters

#endif  // SRC_READ_IN_PARAMETERS_H_
