#ifndef H_michi2
#define H_michi2

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>    // std::transform
//#include <math.h>
#include <cmath>        // std::nan, std::isnan
#include <pthread.h>
#include <unistd.h>
#include <vector>
#include <string>
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw, std::setprecision
#include <sstream>      //std::istringstream
#include <fstream>      // std::ifstream, std::ofstream
#include "michi2_DataClass.h"
#include "michi2_MinPack.h"
#include "michi2_Constraint.h"
//#include "spline.cpp"
//#include "integrate.cpp"
//#include "currentdatetime.cpp"

using namespace std;



extern const char *InfoRedshift;

extern int NumbParallel; // number of parallel subprocesses




//extern const std::string currentDateTime();
//extern vector<double> spline(vector<double> &x, vector<double> &y, vector<double> &output_x);
//extern double integrate(vector<double> &x, vector<double> &y, vector<double> &xrange);
//extern double integrate_LIR(vector<double> &x, vector<double> &y, vector<double> &xrange);



extern std::vector<michi2Constraint *> Constraints;

typedef struct {std::vector<double> X; std::vector<double> Y; std::string Name;} FilterCurveXY;

extern std::vector<FilterCurveXY *> FilterCurves;



typedef struct {std::vector<double> j0; std::vector<double> f0; std::vector<double> df; std::vector<double> f1;} MatchedObsLibStruct;

MatchedObsLibStruct michi2MatchObs(michi2DataClass *DCOBS, michi2DataClass *DCLIB, int debug = 0);



double michi2GetChiSquare(std::vector<double> f1, std::vector<double> f0, std::vector<double> df, double *a1);

double michi2GetReducedChiSquare(std::vector<double> f1, std::vector<double> f0, std::vector<double> df, double *a1);

double michi2VecMean(std::vector<double> vec);

void mnchi2(std::vector<std::string> InputObsList, std::vector<std::string> InputLibList, std::vector<std::string> OutputTableList, std::vector<std::string> InputFilterCurveList);

void *mnchi2parallel(void *params);









#endif
