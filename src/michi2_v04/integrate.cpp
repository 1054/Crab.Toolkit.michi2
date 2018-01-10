#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;


/*
 
 calculate the integration of a curve within xrange
 here the integration is not simply addition, but is \int f_nu dnu
 
 */



double integrate(vector<double> &x, vector<double> &y, vector<double> &xrange)
{
    
    // x must be in increasing order
    
    long i=0, j=0, k=0;
    double output_integral = 0.0;
    double output_width = 0.0; // calc the width of added x data in x axis unit
    long output_count = 0; // count how many x data added
    for(k=0; k<xrange.size(); k+=2) {
        double xlower = (xrange[k+1]<xrange[k+0]) ? xrange[k+1]: xrange[k+0];
        double xupper = (xrange[k+1]<xrange[k+0]) ? xrange[k+0]: xrange[k+1];
        for(j=0; j<x.size(); j++) {
            if(x[j]<xlower) {continue;}
            if(x[j]>xupper) {break;}
            output_integral += y[j];
            output_count ++;
            if(j==0) {
                output_width += x[j+1]-x[j];
            } else {
                output_width += x[j]-x[j-1];
            }
        }
    }
    
    return output_integral;
}



double integrate_LIR(vector<double> &x, vector<double> &y, vector<double> &xrange)
{
    
    // x must be in increasing order
    
    long i=0, j=0, k=0;
    double output_integral = 0.0;
    double output_width = 0.0; // calc the width of added x data in x axis unit
    long output_count = 0; // count how many x data added
    for(k=0; k<xrange.size(); k+=2) {
        double xlower = (xrange[k+1]<xrange[k+0]) ? xrange[k+1]: xrange[k+0];
        double xupper = (xrange[k+1]<xrange[k+0]) ? xrange[k+0]: xrange[k+1];
        for(j=0; j<x.size(); j++) {
            if(x[j]<xlower) {continue;}
            if(x[j]>xupper) {break;}
            double LIR_freq_width = 0.0; // calc the width of current x data during the loop
            if(j==0) {
                LIR_freq_width = fabs(2.99792458e5/x[j+1]-2.99792458e5/x[j]); // assuming x unit is um, this width has a unit of GHz.
                output_width += x[j+1]-x[j];
            } else {
                LIR_freq_width = fabs(2.99792458e5/x[j]-2.99792458e5/x[j-1]); // assuming x unit is um, this width has a unit of GHz.
                output_width += x[j]-x[j-1];
            }
            output_integral += y[j] * LIR_freq_width; // assuming y unit is mJy, this integral has a unit of mJy GHz.
            output_count ++;
        }
    }
    
    return output_integral;
}











