/*
 
 This code aims at testing spline.cpp
 
 to compile and run the test:
    clang++ -I3rd/epsplot -Wno-format-security -std=c++11 test_spline.cpp -o test_spline
    ./test_spline
 
 */


#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include "spline.cpp"
#include "3rd/epsplot/EPSPlot.h"

using namespace std;

extern vector<double> spline(vector<double> &x, vector<double> &y, vector<double> &output_x);


/*
 
class EPSPlotRect {
public:
    double x1, x2, y1, y2, w, h, xscale, yscale;
    double xrange[2], yrange[2];
    EPSPlotRect(double xrange_low, double xrange_high, double yrange_low, double yrange_high, double x_left = 35, double y_bottom = 35, double width = 500, double height = 500) {
        this->xrange[0]=xrange_low;
        this->yrange[0]=yrange_low;
        this->xrange[1]=xrange_high;
        this->yrange[1]=yrange_high;
        this->x1=x_left;
        this->y1=y_bottom;
        this->w=width;
        this->h=height;
        this->x2=this->x1+this->w;
        this->y2=this->y1+this->h;
        this->xscale = this->w / (xrange[1]-xrange[0]);
        this->yscale = this->h / (yrange[1]-yrange[0]);
    };
    double scaled_x(double x) {return (x-this->xrange[0])*this->xscale+this->x1;};
    double scaled_y(double y) {return (y-this->yrange[0])*this->yscale+this->y1;};
};

class EPSPlotSymbol {
public:
    std::vector<float> x, y;
    EPSPlotSymbol(EPSPlotRect *rect, double xvalue, double yvalue, double symsize = 1.0, std::string shape = "square") {
        if(shape.compare("square")>=0) {
            double xpoint = rect->scaled_x(xvalue);
            double ypoint = rect->scaled_y(yvalue);
            double xsize = rect->w/60.0 * symsize;
            double ysize = rect->h/60.0 * symsize;
            x.push_back(xpoint-xsize/2.0); y.push_back(ypoint);
            x.push_back(xpoint-xsize/2.0); y.push_back(ypoint-ysize/2.0);
            x.push_back(xpoint+xsize/2.0); y.push_back(ypoint-ysize/2.0);
            x.push_back(xpoint+xsize/2.0); y.push_back(ypoint+ysize/2.0);
            x.push_back(xpoint-xsize/2.0); y.push_back(ypoint+ysize/2.0);
            x.push_back(xpoint-xsize/2.0); y.push_back(ypoint);
            //<DEBUG>std::cout << xpoint+xsize/2.0 << " " << ypoint+ysize/2.0 << std::endl;
            //<DEBUG>std::cout << xpoint-xsize/2.0 << " " << ypoint-ysize/2.0 << std::endl;
        }
    };
    float *xpoints() {return (float *)(&this->x[0]);} // https://stackoverflow.com/questions/4763590/how-can-i-convert-a-stdvectorfloat-to-a-float-array
    float *ypoints() {return (float *)(&this->y[0]);} // https://stackoverflow.com/questions/4763590/how-can-i-convert-a-stdvectorfloat-to-a-float-array
    int npoints() {return this->x.size();}
};

 */

















int main(int argc, char **argv) {
    
    vector<double> x; x.clear();
    vector<double> y; y.clear();
    vector<double> output_x;
    vector<double> output_y;
    double start = -2.6;
    double end = 2.6;
    double step = 0.05;
    
    
    double val = start;
    while (val < end) {
        x.push_back(val);
        //x.push_back(pow(10,val));
        val += step;
    }
    
    int i=0;
    for(i=0; i<x.size(); i++) {
        y.push_back(2.0*sin(x[i])-0.5*x[i]);
        //y.push_back(2.0*sin(x[i]));
    }
    
    output_x.push_back(0.5);
    output_x.push_back(0.75);
    output_x.push_back(1.25);
    output_x.push_back(3.25);
    
    output_y = spline(x,y,output_x);
    
    
    for(i=0; i<output_x.size(); i++) {
        std::cout << output_x[i] << " " << output_y[i] << std::endl;
    }

    
    EPSPlot plot("testpage.eps");
    plot.set_range(-3.5,3.5, -9.0,9.0);
    plot.plot_line(x,y);
    plot.set_color("red");
    plot.plot_symbol(output_x, output_y, "open square");
    
    
    
    
    return 0;

}












