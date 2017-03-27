#ifndef SRC_EMOJI_H_
#define SRC_EMOJI_H_

#include <cstring>

//! this is a class that includes ascii-emoji strings.

//! I took the text from https://github.com/dysfunc/ascii-emoji
//! and http://japaneseemoticons.me
class emoji {
 public:
    emoji() {};
    ~emoji() {};

    string excited(int id) {
        if (id == 1) return("☜(⌒▽⌒)☞");
        else return("Ｏ(≧▽≦)Ｏ");
    }

    string success() {return("(*•̀ᴗ•́*)و ̑̑");}
    string angry() {return("(╬ ಠ益ಠ)");}
    string very_angry() {return("┻━┻ ︵ヽ(`Д´)ﾉ︵┻━┻");}
    string innocent() {return("ʘ‿ʘ");}
    string cute() {return("(｡◕‿◕｡)");}
    string surprise() {return("（　ﾟДﾟ）");}
    string meh() {return("¯\\(°_o)/¯");}
    string running() {return("ε=ε=ε=┌(;*´Д`)ﾉ");}
    string happy() {return("ヽ(´▽`)/");}
    string disagree() {return("٩◔̯◔۶");}
    string discombobulated() {return("⊙﹏⊙");}
    string sad() {return("¯\\_(⊙︿⊙)_/¯");}
    string confused() {return("¿ⓧ_ⓧﮌ");}
    string confused_scratch() {return("(⊙.☉)7");}
    string dear_god_why() {return("щ（ﾟДﾟщ）");}
    string trolling() {return("༼∵༽ ༼⍨༽ ༼⍢༽ ༼⍤༽");}
    string fuck_it() {return("t(-_-t)");}
    string dancing() {return("┌(ㆆ㉨ㆆ)ʃ");}
    string argry_birds() {return("( ఠൠఠ )ﾉ");}
    string not_support() {return("乁( ◔ ౪◔)「      ┑(￣Д ￣)┍");}
    string pointing() {return("(☞ﾟヮﾟ)☞");}
    string zombie() {return("[¬º-°]¬");}
    string dislike() {return("(Ծ‸ Ծ)");}
    string cry_troll() {return("༼ ༎ຶ ෴ ༎ຶ༽");}
    string cat() {return("((ΦωΦ))");}
    
};


#endif  // SRC_EMOJI_H_
