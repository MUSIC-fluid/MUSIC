#ifndef SRC_READ_IN_PARAMETERS_H_
#define SRC_READ_IN_PARAMETERS_H_

#include <cstring>

#include "./data.h"
#include "./util.h"
#include "./emoji.h"
#include "./pretty_ostream.h"

//! This class handles read in parameters
class ReadInParameters {
 private:
    Util *util;
    emoji *emoji_face;
    string input_file;
    pretty_ostream music_message;

 public:
     ReadInParameters();
     ~ReadInParameters();

    void read_in_parameters(InitData *parameter_list, std::string input_file);
    void check_parameters(InitData *parameter_list);

};

#endif  // SRC_READ_IN_PARAMETERS_H_
