#ifndef EPSPLOT_H
#define EPSPLOT_H

#include <stdio.h>
#include "EPSPlotCore.h"
#include "EPSPlotRect.h"
#include "EPSPlotSymbol.h"
#include "EPSPlotColor.h"

using namespace std;

class EPSPlot {
private:
    EPSPlotCore *m_PlotCore;
    EPSPlotRect *m_PlotRect;
    EPSPlotColor *m_PlotColor;
    bool has_plot_axis_x1;
    bool has_plot_axis_x2;
    bool has_plot_axis_y1;
    bool has_plot_axis_y2;
    bool has_plot_frame_box;
public:
    EPSPlot(std::string FileName) {
        m_PlotRect = new EPSPlotRect(-1.0,1.0,-1.0,1.0);
        m_PlotCore = new EPSPlotCore(FileName.c_str(), m_PlotRect->x1, m_PlotRect->y1, m_PlotRect->x2, m_PlotRect->y2);
        m_PlotColor = new EPSPlotColor();
        m_PlotCore->setLineCap(0);
        m_PlotCore->setLineJoin(1);
        m_PlotCore->setLineWidth(2.0f);
        has_plot_axis_x1 = false;
        has_plot_axis_x2 = false;
        has_plot_axis_y1 = false;
        has_plot_axis_y2 = false;
        has_plot_frame_box = false;
    };
    ~EPSPlot() {
        if(!has_plot_frame_box) {plot_frame_box();}
        m_PlotCore->~EPSPlotCore(); //<TODO>
    };
    void set_range(double xrange_low, double xrange_high, double yrange_low, double yrange_high) {
        if(m_PlotRect) {
            m_PlotRect->set_range(xrange_low, xrange_high, yrange_low, yrange_high);
        }
    };
    void set_color(std::string color_str) {
        if(m_PlotRect) {
            m_PlotColor->set_color(color_str);
            m_PlotCore->setRGBColor(m_PlotColor->R,m_PlotColor->G,m_PlotColor->B);
        }
    };
    void plot_line(std::vector<double> x, std::vector<double> y) {
        int i=0;
        if(m_PlotRect) {
            for(i=1; i<x.size(); i++) {
                m_PlotCore->drawLine(m_PlotRect->scaled_x(x[i-1]), m_PlotRect->scaled_y(y[i-1]), m_PlotRect->scaled_x(x[i]), m_PlotRect->scaled_y(y[i]));
            }
        }
    };
    void plot_symbol(std::vector<double> x, std::vector<double> y, std::string shape = "open square") {
        int i=0, j=0;
        if(m_PlotRect) {
            for(i=0; i<x.size(); i++) {
                EPSPlotSymbol symbol(m_PlotRect, x[i], y[i], shape);
                for(j=0; j<symbol.nsep(); j++) {
                    m_PlotCore->drawLines(symbol.xpoints(j), symbol.ypoints(j), symbol.npoints(j));
                }
            }
        }
    };
    void plot_frame_box() {
        if(m_PlotRect) {
            m_PlotCore->setLineWidth(2.0f);
            m_PlotCore->setRGBColor(0,0,0);
            m_PlotCore->drawBox(m_PlotRect->x1, m_PlotRect->y1, m_PlotRect->x2, m_PlotRect->y2);
        }
    };
    void plot_axis_x1() {
        if(m_PlotRect) {
            if(m_PlotRect->xrange[1]!=m_PlotRect->xrange[0]) {
                // determine starting major tick
                double tick_min = m_PlotRect->xrange[0];
                double tick_max = m_PlotRect->xrange[1];
                double tick_maxabs = std::max(std::fabs(tick_min), std::fabs(tick_max));
                long   tick_expint = long(std::log10(tick_maxabs));
                double tick_level = pow(10,(double)tick_expint);
                double tick_interval = (m_PlotRect->xrange[1]-m_PlotRect->xrange[0])/4; //<TODO> assuming 4 major ticks
                //double xtick_majstart = pow(10,double(long(std::log10(std::fabs(m_PlotRect->xrange[0]))))) * (m_PlotRect->xrange[0]/fabs(m_PlotRect->xrange[0])); //pow(10,double(long(xtick_exponent)));
                // draw major tick
                /*
                 <TODO>
                 <TODO>
                 <TODO>
                 <TODO>
                 <TODO>
                 <TODO>
                double xtick_painter = xtick_majstart;
                while (xtick_painter <= m_PlotRect->xrange[1]) {
                    m_PlotCore->setLineWidth(2.0f);
                    m_PlotCore->setRGBColor(0,0,0);
                    m_PlotCore->drawLine(m_PlotRect->scaled_x(xrval), m_PlotRect->scaled_y(yrval));
                    xrval += (m_PlotRect->xrange[1]-m_PlotRect->xrange[0]);
                }
                 */
            }
        }
    };
    void plot_axis_x2() {
        if(m_PlotRect) {
            m_PlotCore->setLineWidth(2.0f);
            m_PlotCore->setRGBColor(0,0,0);
            /*
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             */
        }
    };
    void plot_axis_y1() {
        if(m_PlotRect) {
            m_PlotCore->setLineWidth(2.0f);
            m_PlotCore->setRGBColor(0,0,0);
            /*
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             */
        }
    };
    void plot_axis_y2() {
        if(m_PlotRect) {
            m_PlotCore->setLineWidth(2.0f);
            m_PlotCore->setRGBColor(0,0,0);
            /*
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             <TODO>
             */
        }
    };
};

#endif
