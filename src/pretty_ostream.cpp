// Copyright Chun Shen @ 2017
// This class is inspired by the JetScapeLogger class written by Joern Putschke

#include "pretty_ostream.h"

#include <sys/resource.h>
#include <sys/time.h>

#include <algorithm>
#include <iomanip>

#include "emoji.h"

using std::cout;
using std::endl;
using std::string;

pretty_ostream::pretty_ostream() {}

pretty_ostream::pretty_ostream(const pretty_ostream&) {}

pretty_ostream& pretty_ostream::operator=(const pretty_ostream&) {
    return (*this);
}

pretty_ostream::~pretty_ostream() {}

std::ostringstream& pretty_ostream::threadBuffer() {
    std::lock_guard<std::mutex> lock(buffersMutex_);
    // std::map never invalidates references to existing elements, and only
    // the owning thread uses its buffer, so it can be used after unlocking
    std::unique_ptr<std::ostringstream>& buffer =
        buffers_[std::this_thread::get_id()];
    if (!buffer) buffer.reset(new std::ostringstream);
    return (*buffer);
}

void pretty_ostream::writeLine(const string& line) {
    static std::mutex outputMutex;
    std::lock_guard<std::mutex> lock(outputMutex);
    cout << line << endl;
}

//! This function flushes out message to the screen
void pretty_ostream::flush(string type) {
    std::ostringstream& buffer = threadBuffer();
    const string message = buffer.str();
    buffer.str("");
    buffer.clear();

    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if (type == "info") {
        info(message);
    } else if (type == "warning") {
        warning(message);
    } else if (type == "error") {
        error(message);
    } else if (type == "debug") {
        debug(message);
    }
}

//! This function output information message
void pretty_ostream::info(string message) {
    std::ostringstream line;
    line << emoji::music_note() << " " << get_memory_usage() << " "
         << message;
    writeLine(line.str());
}

//! This function output debug message
void pretty_ostream::debug(string message) {
    std::ostringstream line;
    line << CYAN << emoji::debug() << " " << get_memory_usage() << " "
         << message << RESET;
    writeLine(line.str());
}

//! This function output warning message
void pretty_ostream::warning(string message) {
    std::ostringstream line;
    line << BOLD << YELLOW << emoji::warning() << " " << message << RESET;
    writeLine(line.str());
}

//! This function output error message
void pretty_ostream::error(string message) {
    std::ostringstream line;
    line << BOLD << RED << emoji::error() << " " << message << RESET;
    writeLine(line.str());
}

//! This function returns a string for the memory usage
//! of the current program in MB
string pretty_ostream::get_memory_usage() {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        double memory_usage_in_MB = 0.0;
#ifdef APPLE
        memory_usage_in_MB = usage.ru_maxrss / 1024. / 1024.;  // MB in Apple
#else
        memory_usage_in_MB = usage.ru_maxrss / 1024.;  // MB in linux
#endif
        std::ostringstream memory_usage;
        memory_usage << std::setprecision(4) << memory_usage_in_MB << " MB";
        return (memory_usage.str());
    } else {
        return ("0 MB");
    }
}
