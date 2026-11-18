// Copyright Chun Shen @ 2017
// This class is inspired by the JetScapeLogger class written by Joern Putschke

#define BOLD "\033[1m"      // Bold
#define BLACK "\033[30m"    // Black
#define RED "\033[31m"      // Red
#define GREEN "\033[32m"    // Green
#define YELLOW "\033[33m"   // Yellow
#define BLUE "\033[34m"     // Blue
#define MAGENTA "\033[35m"  // Magenta
#define CYAN "\033[36m"     // Cyan
#define WHITE "\033[37m"    // White
#define RESET "\033[0m"     // reset

#ifndef PRETTY_OSTREAM_H_
#define PRETTY_OSTREAM_H_

#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

//! Thread-safe message logger. Loggers are shared class members that are
//! also used inside OpenMP parallel loops, so every thread builds its
//! message in its own buffer, and each finished line is written to the
//! screen in one piece under a lock.
class pretty_ostream {
  private:
    std::mutex buffersMutex_;
    std::map<std::thread::id, std::unique_ptr<std::ostringstream>> buffers_;

    //! the message buffer of the calling thread
    std::ostringstream& threadBuffer();

    //! write one complete line to std::cout
    static void writeLine(const std::string& line);

  public:
    pretty_ostream();
    // a copy starts with empty buffers; pending messages are not copied
    pretty_ostream(const pretty_ostream&);
    pretty_ostream& operator=(const pretty_ostream&);
    ~pretty_ostream();

    void flush(std::string type);

    //! This function output information message
    void info(std::string message);

    //! This function output debug message
    void debug(std::string message);

    //! This function output warning message
    void warning(std::string message);

    //! This function output error message
    void error(std::string message);

    //! This function returns a string for the memory usage
    //! of the current program in MB
    std::string get_memory_usage();

    //! reload the << operator
    template <typename T>
    pretty_ostream& operator<<(T const& value) {
        threadBuffer() << value;
        return (*this);
    }
};

#endif  // PRETTY_OSTREAM_H_
