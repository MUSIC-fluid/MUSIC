#ifndef SRC_EMOJI_H_
#define SRC_EMOJI_H_

#include <string>

//! this is a class that includes ascii-emoji strings.

//! I took the text from https://github.com/dysfunc/ascii-emoji
//! and http://japaneseemoticons.me
namespace emoji {
    std::string excited(int id);
    std::string success          ();
    std::string angry            ();
    std::string very_angry       ();
    std::string innocent         ();
    std::string cute             ();
    std::string surprise         ();
    std::string meh              ();
    std::string running          ();
    std::string happy            ();
    std::string disagree         ();
    std::string discombobulated  ();
    std::string sad              ();
    std::string confused         ();
    std::string confused_scratch ();
    std::string dear_god_why     ();
    std::string trolling         ();
    std::string fuck_it          ();
    std::string dancing          ();
    std::string argry_birds      ();
    std::string not_support      ();
    std::string pointing         ();
    std::string zombie           ();
    std::string dislike          ();
    std::string cry_troll        ();
    std::string cat              ();
}

#endif  // SRC_EMOJI_H_
