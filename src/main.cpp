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
    
    if (argc > 1)
        input_file = *(argv+1);
    else
        input_file = "";

    MUSIC *music_hydro = new MUSIC(&DATA, input_file);
    int running_mode = music_hydro->get_running_mode();

    if (running_mode == 1 || running_mode == 2) {
        music_hydro->initialize_hydro();
        music_hydro->run_hydro();
    } else if (running_mode == 73) {
        music_hydro->output_transport_coefficients();
    }
    return(0);
}  /* main */

