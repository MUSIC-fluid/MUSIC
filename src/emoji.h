#ifndef SRC_EMOJI_H_
#define SRC_EMOJI_H_

#include <string>

//! this is a class that includes ascii-emoji strings.

//! I took the text from https://github.com/dysfunc/ascii-emoji
//! and http://japaneseemoticons.me
namespace emoji {
    std::string success            ();
    std::string angry              ();
    std::string very_angry         ();
    std::string innocent           ();
    std::string cute               ();
    std::string surprise           ();
    std::string meh                ();
    std::string happy              ();
    std::string disagree           ();
    std::string sad                ();
    std::string smile_with_sunglass();
    std::string confused           ();
    std::string dislike            ();
    std::string cry_troll          ();
    std::string cat                ();
    std::string thumbup            ();
    std::string beer               ();
    std::string beerclinking       ();
    std::string waterwave          ();
    std::string not_supprot        ();
    std::string raise_hand         ();
    std::string clock              ();
    std::string thinking           ();
    std::string stopwatch          ();
    std::string information        ();
    std::string warning            ();
    std::string error              ();
    std::string debug              ();
    std::string music_note         ();
}

#endif  // SRC_EMOJI_H_
