#include "music_logo.h"
#include <stdio.h>
#include <sys/stat.h>
#include<iostream>

namespace MUSIC_LOGO {

//! This function prints out the program logo
void display_logo(int selector) {
    switch (selector) {
        case 0:  // 3D Diagonal
            std::cout << "================================================================"   << std::endl;
            std::cout << "|           ____                                               |"   << std::endl;
            std::cout << "|         ,'  , `.               .--.--.      ,---,  ,----..   |"   << std::endl;
            std::cout << "|      ,-+-,.' _ |         ,--, /  /    '. ,`--.' | /   /   \\  |"  << std::endl;
            std::cout << "|   ,-+-. ;   , ||       ,'_ /||  :  /`. / |   :  :|   :     : |"   << std::endl;
            std::cout << "|  ,--.'|'   |  ;|  .--. |  | :;  |  |--`  :   |  '.   |  ;. / |"   << std::endl;
            std::cout << "| |   |  ,', |  ':,'_ /| :  . ||  :  ;_    |   :  |.   ; /--`  |"   << std::endl;
            std::cout << "| |   | /  | |  |||  ' | |  . . \\  \\    `. '   '  ;;   | ;     |" << std::endl;
            std::cout << "| '   | :  | :  |,|  | ' |  | |  `----.   \\|   |  ||   : |     |"  << std::endl;
            std::cout << "| ;   . |  ; |--' :  | | :  ' ;  __ \\  \\  |'   :  ;.   | '___  |" << std::endl;
            std::cout << "| |   : |  | ,    |  ; ' |  | ' /  /`--'  /|   |  ''   ; : .'| |"   << std::endl;
            std::cout << "| |   : '  |/     :  | : ;  ; |'--'.     / '   :  |'   | '/  : |"   << std::endl;
            std::cout << "| ;   | |`-'      '  :  `--'   \\ `--'---'  ;   |.' |   :    /  |"  << std::endl;
            std::cout << "| |   ;/          :  ,      .-./           '---'    \\   \\ .'   |" << std::endl;
            std::cout << "| '---'            `--`----'                         `---`     |"   << std::endl;
            std::cout << "================================================================"   << std::endl;
            break;
        case 1:  // bloody
            std::cout << "==============================================" << std::endl;
            std::cout << "|  ███▄ ▄███▓ █    ██   ██████  ██▓ ▄████▄   |" << std::endl;
            std::cout << "| ▓██▒▀█▀ ██▒ ██  ▓██▒▒██    ▒ ▓██▒▒██▀ ▀█   |" << std::endl;
            std::cout << "| ▓██    ▓██░▓██  ▒██░░ ▓██▄   ▒██▒▒▓█    ▄  |" << std::endl;
            std::cout << "| ▒██    ▒██ ▓▓█  ░██░  ▒   ██▒░██░▒▓▓▄ ▄██▒ |" << std::endl;
            std::cout << "| ▒██▒   ░██▒▒▒█████▓ ▒██████▒▒░██░▒ ▓███▀ ░ |" << std::endl;
            std::cout << "| ░ ▒░   ░  ░░▒▓▒ ▒ ▒ ▒ ▒▓▒ ▒ ░░▓  ░ ░▒ ▒  ░ |" << std::endl;
            std::cout << "| ░  ░      ░░░▒░ ░ ░ ░ ░▒  ░ ░ ▒ ░  ░  ▒    |" << std::endl;
            std::cout << "| ░      ░    ░░░ ░ ░ ░  ░  ░   ▒ ░░         |" << std::endl;
            std::cout << "|        ░      ░           ░   ░  ░ ░       |" << std::endl;
            std::cout << "|                                  ░         |" << std::endl;
            std::cout << "==============================================" << std::endl;
            break;
        case 2:  // Dancing font
            std::cout << "===================================================="          << std::endl;
            std::cout << "|   __  __     _   _   ____                   ____  |"         << std::endl;
            std::cout << "| U|' \\/ '|uU |\"|u| | / __\"| u      ___    U /\"___| |"     << std::endl;
            std::cout << "| \\| |\\/| |/ \\| |\\| |<\\___ \\/      |_\"_|   \\| | u   |" << std::endl;
            std::cout << "|  | |  | |   | |_| | u___) |       | |     | |/__  |"         << std::endl;
            std::cout << "|  |_|  |_|  <<\\___/  |____/>>    U/| |\\u    \\____| |"      << std::endl;
            std::cout << "| <<,-,,-.  (__) )(    )(  (__).-,_|___|_,-._// \\\\  |"       << std::endl;
            std::cout << "|  (./  \\.)     (__)  (__)      \\_)-' '-(_/(__)(__) |"       << std::endl;
            std::cout << "===================================================="          << std::endl;
            break;
        case 3:  // STAR Wars
            std::cout << "====================================================="    << std::endl;
            std::cout << "| .___  ___.  __    __       _______. __    ______  |"    << std::endl;
            std::cout << "| |   \\/   | |  |  |  |     /       ||  |  /      | |"   << std::endl;
            std::cout << "| |  \\  /  | |  |  |  |    |   (----`|  | |  ,----' |"   << std::endl;
            std::cout << "| |  |\\/|  | |  |  |  |     \\   \\    |  | |  |      |" << std::endl;
            std::cout << "| |  |  |  | |  `--'  | .----)   |   |  | |  `----. |"    << std::endl;
            std::cout << "| |__|  |__|  \\______/  |_______/    |__|  \\______| |"  << std::endl;
            std::cout << "====================================================="    << std::endl;
            break;
    }

}

//! This function prints out code desciprtion and copyright information
void display_code_description_and_copyright() {
    std::cout << "MUSIC - a 3+1D viscous relativistic hydrodynamic code for "
              << "heavy ion collisions" << std::endl;
    std::cout << "Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, "
              << "Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen"
              << std::endl;
}

//! This function prints out the welcome message
void welcome_message() {
    srand (time(NULL));
    display_logo(rand()%4);
    display_code_description_and_copyright();
}

}
