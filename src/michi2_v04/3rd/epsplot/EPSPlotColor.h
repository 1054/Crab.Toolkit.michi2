#ifndef EPSPLOTCOLOR_H
#define EPSPLOTCOLOR_H

#include <stdio.h>
#include <algorithm>
#include <string>
#include <cctype>

using namespace std;

// requires C++11

// Modified by dzliu,
// 2017-09-30

class EPSPlotColor {
public:
    float R, G, B;
    
    EPSPlotColor(std::string color_str = "default") {
        R=0.0; G=0.0; B=0.0; // black
        set_color(color_str);
    };
    
    ~EPSPlotColor() {
    };
    
    void set_color(std::string color_str = "default") {
        // http://www.rapidtables.com/web/color/RGB_Color.htm
        if     (matchWholeStringIC(color_str,"default"  )) {R=0  /255.0; G=0  /255.0; B=0  /255.0;}
        else if(matchWholeStringIC(color_str,"white"    )) {R=255/255.0; G=255/255.0; B=255/255.0;}
        else if(matchWholeStringIC(color_str,"Red"      )) {R=255/255.0; B=0  /255.0; G=0  /255.0;}
        else if(matchWholeStringIC(color_str,"Lime"     )) {R=0  /255.0; B=255/255.0; G=0  /255.0;}
        else if(matchWholeStringIC(color_str,"Blue"     )) {R=0  /255.0; B=0  /255.0; G=255/255.0;}
        else if(matchWholeStringIC(color_str,"Yellow"   )) {R=255/255.0; B=255/255.0; G=0  /255.0;}
        else if(matchWholeStringIC(color_str,"Cyan"     )) {R=0  /255.0; B=255/255.0; G=255/255.0;}
        else if(matchWholeStringIC(color_str,"Aqua"     )) {R=0  /255.0; B=255/255.0; G=255/255.0;}
        else if(matchWholeStringIC(color_str,"Magenta"  )) {R=255/255.0; B=0  /255.0; G=255/255.0;}
        else if(matchWholeStringIC(color_str,"Fuchsia"  )) {R=255/255.0; B=0  /255.0; G=255/255.0;}
        else if(matchWholeStringIC(color_str,"Silver"   )) {R=192/255.0; B=192/255.0; G=192/255.0;}
        else if(matchWholeStringIC(color_str,"Gray"     )) {R=128/255.0; B=128/255.0; G=128/255.0;}
        else if(matchWholeStringIC(color_str,"Maroon"   )) {R=128/255.0; B=0  /255.0; G=0  /255.0;}
        else if(matchWholeStringIC(color_str,"Olive"    )) {R=128/255.0; B=128/255.0; G=0  /255.0;}
        else if(matchWholeStringIC(color_str,"Green"    )) {R=0  /255.0; B=128/255.0; G=0  /255.0;}
        else if(matchWholeStringIC(color_str,"Purple"   )) {R=128/255.0; B=0  /255.0; G=128/255.0;}
        else if(matchWholeStringIC(color_str,"Teal"     )) {R=0  /255.0; B=128/255.0; G=128/255.0;}
        else if(matchWholeStringIC(color_str,"Navy"     )) {R=0  /255.0; B=0  /255.0; G=128/255.0;}
    };
    
    
    
    
    
    
    
    
    
    // String compare - Try to find in the Haystack the Needle - ignore case
    // -- https://stackoverflow.com/questions/3152241/case-insensitive-stdstring-find
    bool matchStringIC(const std::string & strHaystack, const std::string & strNeedle) {
        auto it = std::search(
                              strHaystack.begin(), strHaystack.end(),
                              strNeedle.begin(),   strNeedle.end(),
                              [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); }
                              );
        return (it != strHaystack.end() );
    };
    
    bool matchWholeStringIC(const std::string & strHaystack, const std::string & strNeedle) {
        auto it = std::search(
                              strHaystack.begin(), strHaystack.end(),
                              strNeedle.begin(),   strNeedle.end(),
                              [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); }
                              );
        return (it != strHaystack.end() );
    };

};

#endif
