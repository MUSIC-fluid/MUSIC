// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016
#include <stdio.h>
#include <sys/stat.h>

#include "./music.h"

using namespace std;

// main program
int main(int argc, char *argv[]) {
    string input_file;
    InitData DATA;
    // you have the option to give a second command line option,
    // which is an integer to be added to the random seed from the current time
    // because on the cluster it will happen that two jobs start at exactly
    // the same time, this makes sure that they dont run with excatly the same
    // seed
    string sseed;
    int seed = 0;
    if (argc > 2) {
        sseed = argv[2];
        seed = atoi(sseed.c_str());
    }
    seed *= 10000;
    DATA.seed = seed;

    if (argc > 1)
        input_file = *(argv+1);
    else
        input_file = "";

    MUSIC *music_ptr = new MUSIC(&DATA, input_file);

    if (DATA.mode == 1 || DATA.mode == 2) {
        music_ptr->run_hydro();
    }
}  /* main */

