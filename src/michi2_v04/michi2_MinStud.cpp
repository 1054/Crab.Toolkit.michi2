#ifndef H_michi2MinStud
#define H_michi2MinStud
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <numeric>
#include <cmath>
#include "CrabTableReadColumn.cpp"
#include "CrabStringReadColumn.cpp"
#include "CrabTableReadInfo.cpp"
#include "CrabTableGetLineCount.cpp"
#include "cminpack/minpack.h"
#include "cminpack/lmdif1_.c"
#include "cminpack/lmdif_.c"
#include "cminpack/lmpar_.c"
#include "cminpack/dpmpar_.c"
#include "cminpack/enorm_.c"
#include "cminpack/fdjac2_.c"
#include "cminpack/qrfac_.c"
#include "cminpack/qrsolv_.c"

using namespace std;


std::vector<double> michi2MinStud_fOBS;
std::vector<double> michi2MinStud_eOBS;
std::vector<std::vector<double> > michi2MinStud_fLIB;
long michi2MinStud_ncount = 0;


void michi2MinStud_func(const int *m, const int *n, const double *x, double *fvec, int *iflag);

void michi2MinStud_func(const int *m, const int *n, const double *x, double *fvec, int *iflag)
{
    /* calculate the functions at x and return the values in fvec[0] through fvec[m-1] -- from cminStud manual html */
    // m is NVAR[0], in LVG, it's the mol line count, in SED, it's the band count.
    double chi2sum = 0.0;
    for (int iim=0; iim<(*m); iim++) {
        double fsum = 0.0;
        double acoe[*n];
        // for (int iin=0; iin<(*n); iin++) { if(x[iin]<0.0) acoe[iin]=0.0; else acoe[iin]=x[iin]; } // <TODO> how to prevent coeff a to be negative?
        // for (int iin=0; iin<(*n); iin++) { fsum += acoe[iin] * michi2MinStud_fLIB[iin][iim]; } // x is the coefficiency a
        for (int iin=0; iin<(*n); iin++) { fsum += x[iin] * michi2MinStud_fLIB[iin][iim]; } // x is the coefficiency a
        fvec[iim] = fsum - michi2MinStud_fOBS[iim]; // chi
        if(!michi2MinStud_eOBS.empty()) { fvec[iim] = fvec[iim]/michi2MinStud_eOBS[iim]; } // chi weighted by obs err
        fvec[iim] = fvec[iim] * fvec[iim]; // chi-square DO NOT SQUARE THE fvec BY OURSELVES! <Corrected><20140822><DzLIU>
        chi2sum += fvec[iim]*fvec[iim]; // fvec is not chi-square! fvec is the list of m non-linear functions in v variables. lmdif_ will calculate the sum of the square of items in fvec, thus we will not square the fvec by ourselves!
    }
    // std::cout << "michi2MinStud::michi2MinStud_func() ncount=" << michi2MinStud_ncount << " a1=" << x[0] << " chi2=" << chi2sum << std::endl;
    michi2MinStud_ncount++;
}










/* Struct Min Stud */
/*
   based on cminStud by 
 */
class michi2MinStud {
public:
    std::vector<double> fOBS;
    std::vector<double> eOBS;
    std::vector<std::vector<double> > fLIB;
    std::vector<double> aCOE;
    std::vector<double> chi2;
    michi2MinStud(std::vector<std::vector<double> > Input_fLIB, std::vector<double> Input_fOBS, std::vector<double> Input_eOBS);
    void func(const int *m, const int *n, const double *x, double *fvec, int *iflag);
    void init(std::vector<std::vector<double> > Input_fLIB, std::vector<double> Input_fOBS, std::vector<double> Input_eOBS);
    void fit();
    double mean(std::vector<double> data);
};










michi2MinStud::michi2MinStud(std::vector<std::vector<double> > Input_fLIB, std::vector<double> Input_fOBS, std::vector<double> Input_eOBS)
{
    //
    init(Input_fLIB, Input_fOBS, Input_eOBS);
    // 
    fit();
    //
}


void michi2MinStud::init(std::vector<std::vector<double> > Input_fLIB, std::vector<double> Input_fOBS, std::vector<double> Input_eOBS)
{
    // save data
    this->fLIB = Input_fLIB;
    this->fOBS = Input_fOBS;
    this->eOBS = Input_eOBS;
    //
    michi2MinStud_fLIB = this->fLIB;
    michi2MinStud_fOBS = this->fOBS;
    michi2MinStud_eOBS = this->eOBS;
    // initialize fitting parameter vector
    this->aCOE.clear();
    this->aCOE.resize(this->fLIB.size());
    for(int i=0; i<this->aCOE.size(); i++) {
        double vLIB = 0.0;
        double vOBS = 0.0;
        for(int j=0; j<this->fLIB[i].size(); j++) { vLIB += this->fLIB[i][j]; }
        for(int j=0; j<this->fLIB[i].size(); j++) { vOBS += this->fOBS[j]; }
        this->aCOE[i] = vOBS/vLIB; // <TODO> WILL ALWAYS BE POSITIVE! THIS IS FOR ASTROPHYSICS!
        // this->aCOE[i] = 0.1;
    }
    //
    this->chi2.clear();
    this->chi2.resize(this->fOBS.size());
    //
    michi2MinStud_ncount = 0;
}



void michi2MinStud::fit()
{
    //
    int m = this->fOBS.size();     // NVAR[0]  N_independent_variables
    int n = this->fLIB.size();     // NCOE     N_composit_coefficients
    this->aCOE.resize(n);
    double *x = this->aCOE.data(); // parameters vector
    double *fvec = new double[m];  // chi-function vector
    double *chisq = new double[m]; // chi-square vector // the square of fvec
    double ftol = 1e-13;           // tolerance <TODO>
    int info = 9;                  // output fit information see manual html
    int lwa = m*n + 5*n + m;       // for internal use
    double *wa = new double[lwa];  // for internal use
    int *iwa = new int[n];         // for internal use
    //
    lmdif1_( michi2MinStud_func, &m, &n, x, fvec, &ftol, &info, iwa, wa, &lwa);
    for(int i=0; i<m; i++) { chisq[i] = fvec[i]*fvec[i]; }
    //
    // save to class
    this->chi2.resize(m);
    this->chi2.assign(chisq,chisq+m); //<TODO> need to check whether x and aCOE are identical.
}



double michi2MinStud::mean(std::vector<double> data)
{
    // return the mean value of *data
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / (double)data.size();
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / data.size() - mean * mean);
    return mean;
}





#endif
