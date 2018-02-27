#include <string>     
#include "emoji.h"

//! this is a class that includes ascii-emoji strings.

//! I took the text from https://github.com/dysfunc/ascii-emoji
//! and http://japaneseemoticons.me

namespace emoji {

std::string excited(int id) {
    if (id == 1) return("☜(⌒▽⌒)☞");
    else return("Ｏ(≧▽≦)Ｏ");
}

std::string success           () {return("(*•̀ᴗ•́*)و ̑̑");}
std::string angry             () {return("(╬ ಠ益ಠ)");}
std::string very_angry        () {return("┻━┻ ︵ヽ(`Д´)ﾉ︵┻━┻");}
std::string innocent          () {return("ʘ‿ʘ");}
std::string cute              () {return("(｡◕‿◕｡)");}
std::string surprise          () {return("（　ﾟДﾟ）");}
std::string meh               () {return("¯\\(°_o)/¯");}
std::string running           () {return("ε=ε=ε=┌(;*´Д`)ﾉ");}
std::string happy             () {return("ヽ(´▽`)/");}
std::string disagree          () {return("٩◔̯◔۶");}
std::string discombobulated   () {return("⊙﹏⊙");}
std::string sad               () {return("¯\\_(⊙︿⊙)_/¯");}
std::string confused          () {return("¿ⓧ_ⓧﮌ");}
std::string confused_scratch  () {return("(⊙.☉)7");}
std::string dear_god_why      () {return("щ（ﾟДﾟщ）");}
std::string trolling          () {return("༼∵༽ ༼⍨༽ ༼⍢༽ ༼⍤༽");}
std::string fuck_it           () {return("t(-_-t)");}
std::string dancing           () {return("┌(ㆆ㉨ㆆ)ʃ");}
std::string argry_birds       () {return("( ఠൠఠ )ﾉ");}
std::string not_support       () {return("乁( ◔ ౪◔)「      ┑(￣Д ￣)┍");}
std::string pointing          () {return("(☞ﾟヮﾟ)☞");}
std::string zombie            () {return("[¬º-°]¬");}
std::string dislike           () {return("(Ծ‸ Ծ)");}
std::string cry_troll         () {return("༼ ༎ຶ ෴ ༎ຶ༽");}
std::string cat               () {return("((ΦωΦ))");}    

}