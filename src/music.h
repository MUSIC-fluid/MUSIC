// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#ifndef SRC_MUSIC_H_
#define SRC_MUSIC_H_

#include "./util.h"
#include "./grid.h"
#include "./data.h"
#include "./init.h"
#include "./eos.h"
#include "./evolve.h"

class MUSIC {
    // this is wrapper class for MUSIC so that it can be used as a external
    // library for integrated framework, such as JETSCAPE
 private:
     int mode;            // records running mode

     InitData *DATA;

     Util *util;
     EOS *eos;

     Grid ***arena;

     Init *init;
     Evolve *evolve;

 public:
     MUSIC(InitData *DATA_in, string input_file);
     ~MUSIC();

     int run_hydro();
     void ReadInData3(string file);
};

#endif  // SRC_MUSIC_H_
