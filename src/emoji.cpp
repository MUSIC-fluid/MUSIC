#include <string>     
#include "emoji.h"

//! this is a class that includes ascii-emoji strings.

//! I took the text from https://github.com/dysfunc/ascii-emoji
//! and http://japaneseemoticons.me

namespace emoji {

std::string success            () {return("\xF0\x9F\x98\x81");}
std::string angry              () {return("\xF0\x9F\x98\xA4");}
std::string very_angry         () {return("\xF0\x9F\x98\xA1");}
std::string innocent           () {return("\xF0\x9F\x98\x8B");}
std::string cute               () {return("\xF0\x9F\x98\x9D");}
std::string surprise           () {return("\xF0\x9F\x98\xB1");}
std::string meh                () {return("\xF0\x9F\x98\x85");}
std::string happy              () {return("\xF0\x9F\x98\x9C");}
std::string smile_with_sunglass() {return("\xF0\x9F\x98\x8E");}
std::string disagree           () {return("\xF0\x9F\x91\x8E");}
std::string sad                () {return("\xF0\x9F\x98\x9E");}
std::string confused           () {return("\xF0\x9F\x98\xB5");}
std::string dislike            () {return("\xF0\x9F\x98\xB0");}
std::string cry_troll          () {return("\xF0\x9F\x98\xAD");}
std::string cat                () {return("\xF0\x9F\x90\x88");}
std::string thumbup            () {return("\xF0\x9F\x91\x8D");}
std::string beer               () {return("\xF0\x9F\x8D\xBA");}
std::string beerclinking       () {return("\U0001F37B");}
std::string waterwave          () {return("\xF0\x9F\x8C\x8A");}
std::string not_supprot        () {return("\xF0\x9F\x99\x85");}
std::string raise_hand         () {return("\xE2\x9C\x8B");}
std::string clock              () {return("\xF0\x9F\x95\x91");}
std::string thinking           () {return("\U0001F914");}
std::string stopwatch          () {return("\U000023F1");}
std::string information        () {return("\U00002139");}
std::string warning            () {return("\U000026A0");}
std::string error              () {return("\U0001F6AB");}
std::string debug              () {return("\U0001F50D");}
std::string music_note         () {return("\U0001F3B5");}

}
