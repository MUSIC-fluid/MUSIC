// Copyright Chun Shen @ 2017
// This class is inspired by the JetScapeLogger class written by Joern Putschke

#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <algorithm>

#include "pretty_ostream.h"

using namespace std;

pretty_ostream::pretty_ostream() {}

pretty_ostream::~pretty_ostream() {}


//! This function flushes out message to the screen
void pretty_ostream::flush(string type) {
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if (type == "info") {
        info(message_stream.str());
    } else if (type == "warning") {
        warning(message_stream.str());
    } else if (type == "error") {
        error(message_stream.str());
    } else if (type == "debug") {
        debug(message_stream.str());
    }
    message_stream.str("");
    message_stream.clear();
}

//! This function output information message
void pretty_ostream::info(string message) {
    cout << "[Info] " << get_memory_usage() << " " << message << endl;
}


//! This function output debug message
void pretty_ostream::debug(string message) {
    cout << CYAN << "[debug] " << get_memory_usage() << " "
         << message << RESET << endl;
}


//! This function output warning message
void pretty_ostream::warning(string message) {
    cout << BOLD << YELLOW << "[Warning] " << message << RESET << endl;
}


//! This function output error message
void pretty_ostream::error(string message) {
    cout << BOLD << RED << "[Error] " << message << RESET << endl;
}

//! This function returns a string for the memory usage
//! of the current program in MB
string pretty_ostream::get_memory_usage() {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        double memory_usage_in_MB = 0.0;
#ifdef APPLE
        memory_usage_in_MB = usage.ru_maxrss/1024./1024.;  // MB in Apple
#else
        memory_usage_in_MB = usage.ru_maxrss/1024.;   // MB in linux
#endif
        ostringstream memory_usage;
        memory_usage << setprecision(4)
                     << memory_usage_in_MB << " MB";
        return(memory_usage.str());
    } else {
        return(0);
    }
}
