// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016

#include <stdio.h>
#include <sys/stat.h>

#include <iostream>
#include <string>

#include "music.h"
#include "music_logo.h"
#include "version.h"

// main program
int main(int argc, char *argv[]) {
    std::string input_file;
    InitData DATA __attribute__((aligned(64)));

    if (argc > 1)
        input_file = *(argv + 1);
    else
        input_file = "";

    MUSIC_LOGO::welcome_message();
    std::cout << "Version: (git branch:" << GIT_BRANCH
              << ") Commit:" << GIT_COMMIT_HASH << std::endl;

    int ireRunHydroCount = 0;
    bool bReRunHydro = false;
    do {
        MUSIC music_hydro(input_file);
        music_hydro.setReRunCount(ireRunHydroCount);
        music_hydro.setReRunHydro(false);

        int running_mode = music_hydro.get_running_mode();
        if (running_mode == 1 || running_mode == 2) {
            music_hydro.initialize_hydro();
            music_hydro.run_hydro();
            bReRunHydro = music_hydro.getReRunHydro();
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
        ireRunHydroCount++;
    } while (bReRunHydro && ireRunHydroCount < 10);
    return (0);
} /* main */
