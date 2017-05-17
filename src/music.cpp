// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include "./music.h"
#include "./dissipative.h"

using namespace std;

MUSIC::MUSIC(InitData *DATA_in, string input_file) {
    welcome_message();
    DATA = DATA_in;
    reader.read_in_parameters(DATA, input_file);
    mode = DATA->mode;
    eos = new EOS(DATA);
    util = new Util();
    if (DATA->Initial_profile == 12 || DATA->Initial_profile == 13) {
        hydro_source_ptr = new hydro_source(DATA_in);
    } else if (DATA->Initial_profile == 30) {
        hydro_source_ptr = new hydro_source(DATA_in);
    }
    flag_hydro_run = 0;
    flag_hydro_initialized = 0;
}


MUSIC::~MUSIC() {
    if (flag_hydro_initialized == 1) {
        delete init;
    }
    if (flag_hydro_run == 1) {
        delete evolve;
    }
    delete eos;
    delete util;
    if (DATA->Initial_profile == 12 || DATA->Initial_profile == 13) {
        delete hydro_source_ptr;
    } else if (DATA->Initial_profile == 30) {
        delete hydro_source_ptr;
    }
}


//! This function initialize hydro
int MUSIC::initialize_hydro() {
    // clean all the surface files
    system("rm surface.dat surface?.dat surface??.dat 2> /dev/null");

    init = new Init(eos, DATA, hydro_source_ptr);
    init->InitArena(DATA, &arena);
    flag_hydro_initialized = 1;
    return(0);
}


//! this is a shell function to run hydro
int MUSIC::run_hydro() {
    evolve = new Evolve(eos, DATA, hydro_source_ptr);
    evolve->EvolveIt(DATA, arena);
    flag_hydro_run = 1;
    return(0);
}


//! this is a test function to output the transport coefficients as
//! function of T and mu_B
void MUSIC::output_transport_coefficients() {
    music_message << "output transport coefficients as functions of T and muB";
    music_message.flush("info");
    Diss *temp_dissipative_ptr = new Diss(eos, DATA);
    temp_dissipative_ptr->output_eta_over_s_T_and_muB_dependence(DATA);
    temp_dissipative_ptr->output_eta_over_s_along_const_sovernB(DATA);
    temp_dissipative_ptr->output_kappa_T_and_muB_dependence(DATA);
    temp_dissipative_ptr->output_kappa_along_const_sovernB(DATA);
    delete temp_dissipative_ptr;
}


//! This function prints out the program logo
void MUSIC::display_logo(int selector) {
    switch (selector) {
        case 0:  // 3D Diagonal
            cout << "================================================================" << endl;
            cout << "|           ____                                               |" << endl;
            cout << "|         ,'  , `.               .--.--.      ,---,  ,----..   |" << endl;
            cout << "|      ,-+-,.' _ |         ,--, /  /    '. ,`--.' | /   /   \\  |" << endl;
            cout << "|   ,-+-. ;   , ||       ,'_ /||  :  /`. / |   :  :|   :     : |" << endl;
            cout << "|  ,--.'|'   |  ;|  .--. |  | :;  |  |--`  :   |  '.   |  ;. / |" << endl;
            cout << "| |   |  ,', |  ':,'_ /| :  . ||  :  ;_    |   :  |.   ; /--`  |" << endl;
            cout << "| |   | /  | |  |||  ' | |  . . \\  \\    `. '   '  ;;   | ;     |" << endl;
            cout << "| '   | :  | :  |,|  | ' |  | |  `----.   \\|   |  ||   : |     |" << endl;
            cout << "| ;   . |  ; |--' :  | | :  ' ;  __ \\  \\  |'   :  ;.   | '___  |" << endl;
            cout << "| |   : |  | ,    |  ; ' |  | ' /  /`--'  /|   |  ''   ; : .'| |" << endl;
            cout << "| |   : '  |/     :  | : ;  ; |'--'.     / '   :  |'   | '/  : |" << endl;
            cout << "| ;   | |`-'      '  :  `--'   \\ `--'---'  ;   |.' |   :    /  |" << endl;
            cout << "| |   ;/          :  ,      .-./           '---'    \\   \\ .'   |" << endl;
            cout << "| '---'            `--`----'                         `---`     |" << endl;
            cout << "================================================================" << endl;
            break;
        case 1:  // bloody
            cout << "==============================================" << endl;
            cout << "|  ███▄ ▄███▓ █    ██   ██████  ██▓ ▄████▄   |" << endl;
            cout << "| ▓██▒▀█▀ ██▒ ██  ▓██▒▒██    ▒ ▓██▒▒██▀ ▀█   |" << endl;
            cout << "| ▓██    ▓██░▓██  ▒██░░ ▓██▄   ▒██▒▒▓█    ▄  |" << endl;
            cout << "| ▒██    ▒██ ▓▓█  ░██░  ▒   ██▒░██░▒▓▓▄ ▄██▒ |" << endl;
            cout << "| ▒██▒   ░██▒▒▒█████▓ ▒██████▒▒░██░▒ ▓███▀ ░ |" << endl;
            cout << "| ░ ▒░   ░  ░░▒▓▒ ▒ ▒ ▒ ▒▓▒ ▒ ░░▓  ░ ░▒ ▒  ░ |" << endl;
            cout << "| ░  ░      ░░░▒░ ░ ░ ░ ░▒  ░ ░ ▒ ░  ░  ▒    |" << endl;
            cout << "| ░      ░    ░░░ ░ ░ ░  ░  ░   ▒ ░░         |" << endl;
            cout << "|        ░      ░           ░   ░  ░ ░       |" << endl;
            cout << "|                                  ░         |" << endl;
            cout << "==============================================" << endl;
            break;
        case 2:  // Dancing font
            cout << "====================================================" << endl;
            cout << "|   __  __     _   _   ____                   ____  |" << endl;
            cout << "| U|' \\/ '|uU |\"|u| | / __\"| u      ___    U /\"___| |" << endl;
            cout << "| \\| |\\/| |/ \\| |\\| |<\\___ \\/      |_\"_|   \\| | u   |" << endl;
            cout << "|  | |  | |   | |_| | u___) |       | |     | |/__  |" << endl;
            cout << "|  |_|  |_|  <<\\___/  |____/>>    U/| |\\u    \\____| |" << endl;
            cout << "| <<,-,,-.  (__) )(    )(  (__).-,_|___|_,-._// \\\\  |" << endl;
            cout << "|  (./  \\.)     (__)  (__)      \\_)-' '-(_/(__)(__) |" << endl;
            cout << "====================================================" << endl;
            break;
        case 3:  // STAR Wars
            cout << "=====================================================" << endl;
            cout << "| .___  ___.  __    __       _______. __    ______  |" << endl;
            cout << "| |   \\/   | |  |  |  |     /       ||  |  /      | |" << endl;
            cout << "| |  \\  /  | |  |  |  |    |   (----`|  | |  ,----' |" << endl;
            cout << "| |  |\\/|  | |  |  |  |     \\   \\    |  | |  |      |" << endl;
            cout << "| |  |  |  | |  `--'  | .----)   |   |  | |  `----. |" << endl;
            cout << "| |__|  |__|  \\______/  |_______/    |__|  \\______| |" << endl;
            cout << "=====================================================" << endl;
            break;
    }

}


//! This function prints out code desciprtion and copyright information
void MUSIC::display_code_description_and_copyright() {
    cout << "MUSIC - a 3+1D viscous relativistic hydrodynamic code for "
         << "heavy ion collisions" << endl;
    cout << "Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, "
         << "Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen"
         << endl;
}

//! This function prints out the welcome message
void MUSIC::welcome_message() {
    srand (time(NULL));
    display_logo(rand()%4);
    display_code_description_and_copyright();
}
