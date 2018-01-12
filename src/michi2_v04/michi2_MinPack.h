#ifndef H_michi2MinPack
#define H_michi2MinPack

#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>
#include <cmath>
//#include "CrabTableReadColumn.cpp"
//#include "CrabStringReadColumn.cpp"
//#include "CrabTableReadInfo.cpp"
//#include "CrabTableGetLineCount.cpp"

using namespace std;



/* global variables */

extern std::vector<double> michi2MinPack_fOBS;

extern std::vector<double> michi2MinPack_eOBS;

extern std::vector<double> michi2MinPack_aCOE;

extern std::vector<std::vector<double> > michi2MinPack_fLIB;

extern long michi2MinPack_ncount;

struct michi2MinPack_constraint {
    int to=-1; std::vector<int> from;
    std::vector<double> multiplication;
    double addition=0.0;
}; //<Added><20171001> allow to constraint aCOE of LIB, for example, lock LIB1 normalization coefficient = LIB2 normalization coefficient * 100, this is useful in locking radio SED to IR(8-1000um) via the IR-radio correlation.

extern std::vector<michi2MinPack_constraint *> michi2MinPack_constraints; //<Added><20171001>



/* global function */

void michi2MinPack_func(const int *m, const int *n, const double *x, double *fvec, int *iflag);



/* 
   Struct Min Pack
   based on cminpack
*/

class michi2MinPack {
public:
    std::vector<double> fOBS;
    std::vector<double> eOBS;
    std::vector<std::vector<double> > fLIB;
    std::vector<double> aCOE;
    std::vector<double> chi2;
    michi2MinPack(
                  std::vector<std::vector<double> > Input_fLIB,
                  std::vector<double> Input_fOBS,
                  std::vector<double> Input_eOBS,
                  int Input_debug = 0
                  );
    void init(
              std::vector<std::vector<double> > Input_fLIB,
              std::vector<double> Input_fOBS,
              std::vector<double> Input_eOBS,
              int Input_debug = 0
              );
    void constrain(
                   int toLib,
                   std::vector<int> fromLibs,
                   std::vector<double> multiplication_factors
                   );
    void constrain(
                   int toLib,
                   int fromLib,
                   double multiplication_factor
                   );
    void func(const int *m, const int *n, const double *x, double *fvec, int *iflag);
    void fit(int Input_debug = 0);
    double mean(std::vector<double> data);
};











#endif
