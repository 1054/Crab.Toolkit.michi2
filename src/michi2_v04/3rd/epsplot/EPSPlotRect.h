#ifndef EPSPLOTRECT_H
#define EPSPLOTRECT_H

#include <stdio.h>
#include "EPSPlotCore.h"

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
    void set_xrange(double xrange_low, double xrange_high) {
        this->xrange[0]=xrange_low;
        this->xrange[1]=xrange_high;
        this->xscale = this->w / (xrange[1]-xrange[0]);
    };
    void set_yrange(double yrange_low, double yrange_high) {
        this->yrange[0]=yrange_low;
        this->yrange[1]=yrange_high;
        this->yscale = this->h / (yrange[1]-yrange[0]);
    };
    void set_range(double xrange_low, double xrange_high, double yrange_low, double yrange_high) {
        this->set_xrange(xrange_low, xrange_high);
        this->set_yrange(yrange_low, yrange_high);
    };
    double scaled_x(double x) {return (x-this->xrange[0])*this->xscale+this->x1;};
    double scaled_y(double y) {return (y-this->yrange[0])*this->yscale+this->y1;};
    bool isValid() {return (this->xrange[0]!=this->xrange[1] && this->yrange[0]!=this->yrange[1]);}
};

#endif
