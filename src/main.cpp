// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016

#include <dirent.h>
#include <stdio.h>
#include <cstdio>
#include <sys/stat.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "music.h"
#include "music_logo.h"
#include "version.h"

namespace {

//! Remove the *.dat files in the working directory left by a previous run
//! (some outputs are opened in append mode), but keep every input file
//! named in the parameters, e.g. the default freeze-out list
//! eps_freeze_list_s95p_v1.dat or an initial profile given as a bare name.
void removeOldOutputDatFiles(const InitData &DATA) {
    const std::vector<std::string> inputFiles = {
        DATA.initName,      DATA.initName_TA,
        DATA.initName_TB,   DATA.initName_participants,
        DATA.initName_AMPT, DATA.freeze_list_filename,
    };
    // identify input files by device and inode, so relative paths such as
    // "./file.dat" or symbolic links are matched too
    std::vector<std::pair<dev_t, ino_t>> keep;
    for (const auto &name : inputFiles) {
        struct stat st;
        if (!name.empty() && stat(name.c_str(), &st) == 0) {
            keep.emplace_back(st.st_dev, st.st_ino);
        }
    }

    DIR *dir = opendir(".");
    if (dir == nullptr) return;
    while (struct dirent *entry = readdir(dir)) {
        const std::string name = entry->d_name;
        if (name.size() <= 4 || name.compare(name.size() - 4, 4, ".dat") != 0)
            continue;
        struct stat st;
        if (lstat(name.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) continue;
        bool isInput = false;
        for (const auto &k : keep) {
            if (k.first == st.st_dev && k.second == st.st_ino) {
                isInput = true;
                break;
            }
        }
        if (!isInput) std::remove(name.c_str());
    }
    closedir(dir);
}

}  // namespace

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

        removeOldOutputDatFiles(music_hydro.getParameters());
        int running_mode = music_hydro.get_running_mode();
        int evoStatus = 0;
        if (running_mode == 1 || running_mode == 2) {
            music_hydro.initialize_hydro();
            evoStatus = music_hydro.run_hydro();
            bReRunHydro = music_hydro.getReRunHydro();
        }

        if (evoStatus == 1) {
            if (running_mode == 1 || running_mode == 3 || running_mode == 4
                || running_mode == 13 || running_mode == 14) {
                music_hydro.run_Cooper_Frye();
            }
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
