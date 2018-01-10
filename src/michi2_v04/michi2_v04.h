#ifndef H_michi2
#define H_michi2

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>    // std::transform
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <vector>
#include <string>
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw, std::setprecision
#include <sstream>      //std::istringstream
#include <fstream>      // std::ifstream, std::ofstream
#include <time.h>
#include "CrabTableReadColumn.cpp"
#include "CrabStringReadColumn.cpp"
#include "CrabTableReadInfo.cpp"
#include "michi2_DataClass.cpp"
#include "michi2_MinPack.cpp"
#include "michi2_Constraint.h"
#include "spline.cpp"
#include "integrate.cpp"


using namespace std;



const char *InfoRedshift = "";

int NumbParallel = 2; // number of parallel subprocesses

const std::string currentDateTime() {
    // Get current date/time, format is YYYY-MM-DD.HH:mm:ss
    // http://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X %Z", &tstruct);
    return buf;
};

std::vector<michi2Constraint *> Constraints; // --> ./michi2_v04_201607 -obs test/... terminated by signal SIGSEGV (Address boundary error)

typedef struct {std::vector<double> X; std::vector<double> Y; std::string Name;} FilterCurveXY;

std::vector<FilterCurveXY *> FilterCurves;





extern vector<double> spline(vector<double> &x, vector<double> &y, vector<double> &output_x);

extern double integrate(vector<double> &x, vector<double> &y, vector<double> &xrange);

extern double integrate_LIR(vector<double> &x, vector<double> &y, vector<double> &xrange);





typedef struct {std::vector<double> j0; std::vector<double> f0; std::vector<double> df; std::vector<double> f1;} MatchedObsLibStruct;

MatchedObsLibStruct michi2MatchObs(michi2DataClass *DCOBS, michi2DataClass *DCLIB, int debug = 0);

double michi2GetChiSquare(std::vector<double> f1, std::vector<double> f0, std::vector<double> df, double *a1);

double michi2GetReducedChiSquare(std::vector<double> f1, std::vector<double> f0, std::vector<double> df, double *a1);

double michi2VecMean(std::vector<double> vec);

void mnchi2(std::vector<std::string> InputObsList, std::vector<std::string> InputLibList, std::vector<std::string> OutputTableList, std::vector<std::string> InputFilterCurveList);

void *mnchi2parallel(void *params);









#endif
