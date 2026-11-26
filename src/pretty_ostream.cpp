// Copyright Chun Shen @ 2017
// This class is inspired by the JetScapeLogger class written by Joern Putschke

#include "pretty_ostream.h"

#include <unistd.h>

#if defined(__APPLE__)
#include <sys/resource.h>
// from <libproc.h>, declared here because GCC cannot parse the mach
// headers that <libproc.h> includes
extern "C" int proc_pid_rusage(int pid, int flavor, rusage_info_t* buffer);
#endif

#include <algorithm>
#include <fstream>
#include <iomanip>

#include "emoji.h"

using std::cout;
using std::endl;
using std::string;

pretty_ostream::pretty_ostream() {}

pretty_ostream::pretty_ostream(const string& tag) : tag_(tag) {}

pretty_ostream::pretty_ostream(const pretty_ostream& other)
    : tag_(other.tag_) {}

pretty_ostream& pretty_ostream::operator=(const pretty_ostream& other) {
    tag_ = other.tag_;
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
    line << emoji::music_note() << " " << get_memory_usage() << " ";
    if (!tag_.empty()) line << tag_ << ": ";
    line << message;
    writeLine(line.str());
}

//! This function output debug message
void pretty_ostream::debug(string message) {
    std::ostringstream line;
    line << CYAN << emoji::debug() << " " << get_memory_usage() << " ";
    if (!tag_.empty()) line << tag_ << ": ";
    line << message << RESET;
    writeLine(line.str());
}

//! This function output warning message
void pretty_ostream::warning(string message) {
    std::ostringstream line;
    line << BOLD << YELLOW << emoji::warning() << " ";
    if (!tag_.empty()) line << tag_ << ": ";
    line << message << RESET;
    writeLine(line.str());
}

//! This function output error message
void pretty_ostream::error(string message) {
    std::ostringstream line;
    line << BOLD << RED << emoji::error() << " ";
    if (!tag_.empty()) line << tag_ << ": ";
    line << message << RESET;
    writeLine(line.str());
}

//! This function returns a string for the current memory usage (the
//! resident set size) of the program in MB
string pretty_ostream::get_memory_usage() {
    double bytes = -1.;
#if defined(__APPLE__)
    rusage_info_v2 info;
    if (proc_pid_rusage(
            getpid(), RUSAGE_INFO_V2, reinterpret_cast<rusage_info_t*>(&info))
        == 0)
        bytes = static_cast<double>(info.ri_resident_size);
#else
    // /proc/self/statm: total program size, resident set size (in pages)
    std::ifstream statm("/proc/self/statm");
    long pages = 0, resident = 0;
    if (statm >> pages >> resident)
        bytes = static_cast<double>(resident) * sysconf(_SC_PAGESIZE);
#endif
    if (bytes < 0.) return ("0 MB");
    std::ostringstream memory_usage;
    memory_usage << std::setprecision(4) << bytes / 1024. / 1024. << " MB";
    return (memory_usage.str());
}
