// Copyright Chun Shen @ 2017
// This class is inspired by the JetScapeLogger class written by Joern Putschke

#define BOLD    "\033[1m"       // Bold
#define BLACK   "\033[30m"      // Black
#define RED     "\033[31m"      // Red
#define GREEN   "\033[32m"      // Green
#define YELLOW  "\033[33m"      // Yellow
#define BLUE    "\033[34m"      // Blue
#define MAGENTA "\033[35m"      // Magenta
#define CYAN    "\033[36m"      // Cyan
#define WHITE   "\033[37m"      // White
#define RESET   "\033[0m"       // reset
 
#ifndef PRETTY_OSTREAM_H_
#define PRETTY_OSTREAM_H_

#include <iostream>
#include <sstream>
#include <string>

class pretty_ostream {
 private:
    std::ostringstream message_stream;

 public:
    pretty_ostream();
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
    template <typename T> pretty_ostream& operator<<(T const& value) {
        message_stream << value;
        return(*this);
    }
};

#endif  // PRETTY_OSTREAM_H_
