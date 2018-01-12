/*
 
 This code aims at testing spline.cpp
 
 to compile and run the test:
 clang++ -I3rd/epsplot -Wno-format-security -std=c++11 test_filter_curve.cpp -o test_filter_curve
 ./test_filter_curve
 
 */


#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <sstream>
#include <fstream>
#include <vector>
#include "spline.cpp"
#include "3rd/epsplot/EPSPlot.h"

using namespace std;

extern vector<double> spline(vector<double> &x, vector<double> &y, vector<double> &output_x);

std::vector<std::string> FilterCurves; // Obs Filter Curve List


std::vector<double> DCLIB_X;
std::vector<double> DCLIB_Y;






int main(int argc, char **argv) {
    
    
    
    
    DCLIB_X.push_back(22.5);
    DCLIB_X.push_back(23.5);
    DCLIB_X.push_back(24.5);
    DCLIB_X.push_back(27.5);
    DCLIB_Y.push_back(100);
    DCLIB_Y.push_back(120);
    DCLIB_Y.push_back(144);
    DCLIB_Y.push_back(189);
    
    double f1_li=0.0;
    
    
    FilterCurves.push_back("test_filter_curve_24um.dat");
    
    
    
    
    
    // see whether user has input filter curve ascii files
    // the user should input as many filter curve ascii files as the OBS data point number
    int k = 0;
    if(FilterCurves.size()>k) {
        // <20170930>
        if(FilterCurves[k].compare("none")!=0 || FilterCurves[k].compare("null")!=0 ||
           FilterCurves[k].compare("None")!=0 || FilterCurves[k].compare("Null")!=0) {
            std::ifstream FilterCurveFileStream(FilterCurves[k]);
            if(FilterCurveFileStream.is_open()) {
                std::vector<double> FilterCurveX;
                std::vector<double> FilterCurveY;
                std::string FilterCurveStr;
                while(getline(FilterCurveFileStream,FilterCurveStr)) {
                    if (FilterCurveStr[0]=='#') continue;
                    std::istringstream FilterCurveStrStream(FilterCurveStr);
                    std::string FilterCurveStrStr;
                    std::vector<std::string> FilterCurveStrList { std::istream_iterator<std::string>(FilterCurveStrStream), {} }; // https://www.quora.com/How-do-I-split-a-string-by-space-into-an-array-in-c++
                    if(FilterCurveStrList.size()>=2) {
                        FilterCurveX.push_back((double)atof(FilterCurveStrList[0].c_str()));
                        FilterCurveY.push_back((double)atof(FilterCurveStrList[1].c_str()));
                        cout << "FilterCurveX = " << FilterCurveX[FilterCurveX.size()-1] << ", "
                             << "FilterCurveY = " << FilterCurveY[FilterCurveY.size()-1] << std::endl;
                    }
                }
                if(FilterCurveX.size()>0 && FilterCurveY.size()>0) {
                    // apply filter curve to the LIB data
                    // spline FilterCurveX FilterCurveY LIBX FilterFactorY
                    
                    std::vector<double> FilterCurveY_Matched = spline(FilterCurveX, FilterCurveY, DCLIB_X);
                    double FilterCurveIntegrated1 = 0.0; // \int R(\lambda) f(\lambda) d\lambda
                    double FilterCurveIntegrated2 = 0.0; // \int R(\lambda) d\lambda
                    for(int i1=0; i1<DCLIB_Y.size(); i1++) {
                        if(DCLIB_X[i1]>=FilterCurveX.front() && DCLIB_X[i1]<=FilterCurveX.back()) {
                            FilterCurveIntegrated1 += DCLIB_Y[i1] * FilterCurveY_Matched[i1];
                            FilterCurveIntegrated2 += FilterCurveY_Matched[i1];
                        }
                    }
                    f1_li = FilterCurveIntegrated1 / FilterCurveIntegrated2; // now applied filter curve (aka transmission curve)
                }
            } else {
                cout << std::endl;
                cout << "Error! Failed to read the filter curve file \"" << FilterCurves[k] << "\"!" << std::endl;
                cout << std::endl;
            }
        }
    }
    
    
    int i=0;
    for (i=0; i<DCLIB_X.size(); i++) {
        cout << fixed << showpoint;
        cout << "DCLIB_X = " << std::setprecision(9) << DCLIB_X[i] << ", "
        << "DCLIB_Y = " << std::setprecision(9) << DCLIB_Y[i] << std::endl;
    }
    cout << "f1_li = " << f1_li << std::endl;
    
    return 0;
    
}












