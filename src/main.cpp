// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016

#include <stdio.h>
#include <string>
#include <sys/stat.h>

#include "music.h"
#include "music_logo.h"

// main program
int main(int argc, char *argv[]) {
    std::string input_file;
    InitData DATA __attribute__ ((aligned (64)));

    if (argc > 1)
        input_file = *(argv+1);
    else
        input_file = "";

    MUSIC_LOGO::welcome_message();
    MUSIC music_hydro(input_file);
    int running_mode = music_hydro.get_running_mode();

    if (running_mode == 1 || running_mode == 2) {
        music_hydro.initialize_hydro();
        music_hydro.run_hydro();
    }

    if (running_mode == 1 || running_mode == 3 || running_mode == 4
            || running_mode == 13 || running_mode == 14) {
        music_hydro.run_Cooper_Frye();
    }

    if (running_mode == 71) {
        music_hydro.check_eos();
    }
    if (running_mode == 73) {
        music_hydro.output_transport_coefficients();
    }

    return(0);
}  /* main */

