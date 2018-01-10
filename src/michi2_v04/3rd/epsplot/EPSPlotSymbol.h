#ifndef EPSPLOTSYMBOL_H
#define EPSPLOTSYMBOL_H

#include <stdio.h>
#include "EPSPlotCore.h"
#include "EPSPlotRect.h"

using namespace std;

class EPSPlotSymbol {
public:
    std::vector<float> x, y;
    std::vector<long> start; // the symbol can be described by several separated line groups
    std::vector<long> end; // each line group is marked by start and end.
    
    EPSPlotSymbol(EPSPlotRect *rect, double xvalue, double yvalue, std::string shape = "square", double symsize = 1.0) {
        long i=0;
        if(matchStringIC(shape,"square")) {
            double xpoint = rect->scaled_x(xvalue);
            double ypoint = rect->scaled_y(yvalue);
            double xsize = rect->w/60.0 * symsize;
            double ysize = rect->h/60.0 * symsize;
            start.push_back(i);
            x.push_back(xpoint-xsize/2.0); y.push_back(ypoint); i++;
            x.push_back(xpoint-xsize/2.0); y.push_back(ypoint-ysize/2.0); i++;
            x.push_back(xpoint+xsize/2.0); y.push_back(ypoint-ysize/2.0); i++;
            x.push_back(xpoint+xsize/2.0); y.push_back(ypoint+ysize/2.0); i++;
            x.push_back(xpoint-xsize/2.0); y.push_back(ypoint+ysize/2.0); i++;
            x.push_back(xpoint-xsize/2.0); y.push_back(ypoint); i++;
            x.push_back(NAN); y.push_back(NAN); i++; end.push_back(i); // we use NAN to separate
            //<DEBUG>std::cout << xpoint+xsize/2.0 << " " << ypoint+ysize/2.0 << std::endl;
            //<DEBUG>std::cout << xpoint-xsize/2.0 << " " << ypoint-ysize/2.0 << std::endl;
        }
        //TODO: filled symbol?
    };
    
    float *xpoints(long isep = 0) {
        if(isep<0 || isep>=this->start.size()) {return NULL;}
        return (float *)(&this->x[this->start[isep]]);
    }; // https://stackoverflow.com/questions/4763590/how-can-i-convert-a-stdvectorfloat-to-a-float-array
    
    float *ypoints(long isep = 0) {
        if(isep<0 || isep>=this->start.size()) {return NULL;}
        return (float *)(&this->y[this->start[isep]]);
    }; // https://stackoverflow.com/questions/4763590/how-can-i-convert-a-stdvectorfloat-to-a-float-array
    
    int npoints(long isep = 0) {
        if(isep<0 || isep>=this->start.size()) {return 0;}
        return (this->end[isep] - this->start[isep] - 1); // the end mark is NAN, so the valid point number is (end-start) without +1.
    };
    
    int nsep() {
        return this->start.size();
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
