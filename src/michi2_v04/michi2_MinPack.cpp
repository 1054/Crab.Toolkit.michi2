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


std::vector<double>                 michi2MinPack_fOBS;
std::vector<double>                 michi2MinPack_eOBS;
std::vector<double>                 michi2MinPack_aCOE;
std::vector<std::vector<double> >   michi2MinPack_fLIB;
long                                michi2MinPack_ncount = 0;
struct michi2MinPack_constraint { int to=-1; std::vector<int> from; std::vector<double> multiplication; double addition=0.0; }; //<Added><20171001> allow to constraint aCOE of LIB, for example, lock LIB1 normalization coefficient = LIB2 normalization coefficient * 100, this is useful in locking radio SED to IR(8-1000um) via the IR-radio correlation.
std::vector<michi2MinPack_constraint *> michi2MinPack_constraints; //<Added><20171001>










void michi2MinPack_func(const int *m, const int *n, const double *x, double *fvec, int *iflag);

void michi2MinPack_func(const int *m, const int *n, const double *x, double *fvec, int *iflag)
{
    /* calculate the functions at x and return the values in fvec[0] through fvec[m-1] -- from cminpack manual html */
    // m is NVAR[0], in LVG, it's the mol line count, in SED, it's the band count.
    // n is the LIB number
    int debug = 0;
    double chi2sum = 0.0;
    for (int iim=0; iim<(*m); iim++) {
        double fsum = 0.0; // sum of multiple libs at one band
        double *acoe = (double *)x; // x is the coefficiency a
        for (int iin=0; iin<(*n); iin++) { if(x[iin]<0.0) acoe[iin]=0.0; else acoe[iin]=x[iin]; } // <TODO> how to prevent coeff a to be negative?
        //
        // <added><20171001> apply constraints directly on LIB normalization factors (i.e. aCOE)
        if(michi2MinPack_constraints.size()>0) {
            for(int iicon=0; iicon<michi2MinPack_constraints.size(); iicon++) {
                michi2MinPack_constraint *temp_constraint = michi2MinPack_constraints[iicon];
                if(temp_constraint->to>=1 && temp_constraint->to<=(*n)) {
                    acoe[temp_constraint->to-1] = 0.0;
                    for(int ijcon=0; ijcon<temp_constraint->from.size(); ijcon++) {
                        if(temp_constraint->from[ijcon]>=1 && temp_constraint->from[ijcon]<=(*n)) {
                            acoe[temp_constraint->to-1] += acoe[temp_constraint->from[ijcon]-1] * temp_constraint->multiplication[ijcon];
                        }
                    }
                    acoe[temp_constraint->to-1] = acoe[temp_constraint->to-1] + temp_constraint->addition;
                    //<DEBUG><20171001> acoe[temp_constraint->to-1] = 1.0;
                }
            }
        }
        //
        //for (int iin=0; iin<(*n); iin++) { michi2MinPack_aCOE[iin] = acoe[iin]; } // <TODO> how to prevent coeff a to be negative?
        //
        for (int iin=0; iin<(*n); iin++) {
            if(!std::isnan(michi2MinPack_fLIB[iin][iim])) {
                fsum += acoe[iin] * michi2MinPack_fLIB[iin][iim];
                //<DEBUG><20171001> std::cout << "acoe[" << iin << "] " << acoe[iin] << " * " << michi2MinPack_fLIB[iin][iim] << std::endl; //<DEBUG><20171001>
            }
        }
        //
        //for (int iin=0; iin<(*n); iin++) {
        //    if(!std::isnan(michi2MinPack_fLIB[iin][iim])) {
        //        fsum += x[iin] * michi2MinPack_fLIB[iin][iim];
        //    /* std::cout << x[iin] << "*" << michi2MinPack_fLIB[iin][iim] << std::endl; */
        //    }
        //} // <Corrected><20140922><DzLIU> fLib smaller than fObs!
        //
        fvec[iim] = abs(fsum - michi2MinPack_fOBS[iim]); // chi
        if(!michi2MinPack_eOBS.empty()) { fvec[iim] = fvec[iim]/michi2MinPack_eOBS[iim]; } // chi weighted by obs err
        // fvec[iim] = fvec[iim] * fvec[iim]; // chi-square DO NOT SQUARE THE fvec BY OURSELVES! <Corrected><20140822><DzLIU>
        chi2sum += fvec[iim]*fvec[iim]; // fvec is not chi-square! fvec is the list of m non-linear functions in v variables. lmdif_ will calculate the sum of the square of items in fvec, thus we will not square the fvec by ourselves!
    }
    if(debug) {
        std::cout << "michi2MinPack::michi2MinPack_func() ncount=" << michi2MinPack_ncount << " bands=" << michi2MinPack_fLIB[0].size() << "," << michi2MinPack_fLIB[1].size() << "," << michi2MinPack_fLIB[2].size() <<  " a1=" << x[0] << " a2=" << x[1] << " a3=" << x[2] << " chi2=" << chi2sum << std::endl;
    }
    michi2MinPack_ncount++;
}










/* Struct Min Pack */
/*
   based on cminpack by 
 */
class michi2MinPack {
public:
    std::vector<double> fOBS;
    std::vector<double> eOBS;
    std::vector<std::vector<double> > fLIB;
    std::vector<double> aCOE;
    std::vector<double> chi2;
    michi2MinPack(std::vector<std::vector<double> > Input_fLIB, std::vector<double> Input_fOBS, std::vector<double> Input_eOBS, int Input_debug = 0);
    void func(const int *m, const int *n, const double *x, double *fvec, int *iflag);
    void init(std::vector<std::vector<double> > Input_fLIB, std::vector<double> Input_fOBS, std::vector<double> Input_eOBS, int Input_debug = 0);
    void constrain(int toLib, std::vector<int> fromLibs, std::vector<double> multiplication_factors);
    void constrain(int toLib, int fromLib, double multiplication_factor);
    void fit(int Input_debug = 0);
    double mean(std::vector<double> data);
};










michi2MinPack::michi2MinPack(std::vector<std::vector<double> > Input_fLIB, std::vector<double> Input_fOBS, std::vector<double> Input_eOBS, int Input_debug)
{
    //
    init(Input_fLIB, Input_fOBS, Input_eOBS, Input_debug);
    // 
    //fit(Input_debug);
    //
}


void michi2MinPack::init(std::vector<std::vector<double> > Input_fLIB, std::vector<double> Input_fOBS, std::vector<double> Input_eOBS, int Input_debug)
{
    // save data
    this->fLIB = Input_fLIB;
    this->fOBS = Input_fOBS;
    this->eOBS = Input_eOBS;
    //
    michi2MinPack_fLIB = this->fLIB;
    michi2MinPack_fOBS = this->fOBS;
    michi2MinPack_eOBS = this->eOBS;
    //
    int debug = Input_debug;
    if(debug) { std::cout << "michi2MinPack: initial"; }
    // initialize fitting parameter vector
    this->aCOE.clear();
    this->aCOE.resize(this->fLIB.size());
    michi2MinPack_aCOE.clear();
    michi2MinPack_aCOE.resize(this->fLIB.size());
    for(int i=0; i<this->aCOE.size(); i++) {
//        double vLIB = 0.0;
//        double vOBS = 0.0;
//        for(int j=0; j<this->fLIB[i].size(); j++) { vLIB += this->fLIB[i][j]; }
//        for(int j=0; j<this->fLIB[i].size(); j++) { vOBS += this->fOBS[j]; }
//        this->aCOE[i] = vOBS/vLIB/double(this->aCOE.size());
        double iniCOE = 1.0;
        for(int j=0; j<this->fLIB[i].size(); j++) {
            if(this->fLIB[i][j]>0.0) {
                if(iniCOE==1.0) { iniCOE = (this->fOBS[j]/this->fLIB[i][j]); }
                if((this->fOBS[j]/this->fLIB[i][j])<iniCOE) { iniCOE = (this->fOBS[j]/this->fLIB[i][j]); }
            }
        } /* <Fixed><20141113><DzLIU> */
        this->aCOE[i] = iniCOE;
        michi2MinPack_aCOE[i] = iniCOE;
        // this->aCOE[i] = 0.1;
        // this->aCOE[i] = 100.0;
        if(debug) { std::cout << " a[" << i << "]=" << iniCOE; }
    }
    if(debug) { std::cout << std::endl; }
    //
    this->chi2.clear();
    this->chi2.resize(this->fOBS.size());
    //
    michi2MinPack_ncount = 0;
}



void michi2MinPack::constrain(int toLib, std::vector<int> fromLibs, std::vector<double> multiplication_factors)
{
    // <20171001>
    michi2MinPack_constraint *temp_constraint = new michi2MinPack_constraint();
    temp_constraint->to = toLib;
    temp_constraint->from.resize(fromLibs.size()); for(int k=0; k<fromLibs.size(); k++) {temp_constraint->from[k] = fromLibs[k];}
    temp_constraint->multiplication.resize(multiplication_factors.size()); for(int k=0; k<multiplication_factors.size(); k++) {temp_constraint->multiplication[k] = multiplication_factors[k];}
    michi2MinPack_constraints.push_back(temp_constraint);
}

void michi2MinPack::constrain(int toLib, int fromLib, double multiplication_factor)
{
    // <20171001>
    michi2MinPack_constraint *temp_constraint = new michi2MinPack_constraint();
    temp_constraint->to = toLib;
    temp_constraint->from.resize(1); for(int k=0; k<1; k++) {temp_constraint->from[k] = fromLib;}
    temp_constraint->multiplication.resize(1); for(int k=0; k<1; k++) {temp_constraint->multiplication[k] = multiplication_factor;}
    michi2MinPack_constraints.push_back(temp_constraint);
}



void michi2MinPack::fit(int Input_debug)
{
    //
    int m = this->fOBS.size();     // NVAR[0]  N_independent_variables
    int n = this->fLIB.size();     // NCOE     N_composit_coefficients
    //this->aCOE.resize(n);
    double *x = new double[n]; for(int i=0; i<n; i++) { x[i] = this->aCOE[i]; } // michi2MinPack_aCOE.data(); // this->aCOE.data(); // parameters vector
    double *fvec = new double[m];  // chi-function vector
    double *chisq = new double[m]; // chi-square vector // the square of fvec
    double ftol = 1e-15;           // tolerance <TODO>
    int info = 9;                  // output fit information see manual html
    int lwa = m*n + 16*n + m;      // for internal use
    double *wa = new double[lwa];  // for internal use
    int *iwa = new int[n];         // for internal use
    //
    lmdif1_( michi2MinPack_func, &m, &n, x, fvec, &ftol, &info, iwa, wa, &lwa );
    //
    // get result
    for(int i=0; i<n; i++) { this->aCOE[i] = x[i]; }
    //
    // prevent from negative COE <20171001>
    //for(int i=0; i<n; i++) { this->aCOE[i] = michi2MinPack_aCOE[i]>0 ? michi2MinPack_aCOE[i] : 0.0; }
    //for(int i=0; i<n; i++) { if(this->aCOE[i]<0.0) {this->aCOE[i] = 0.0;}; michi2MinPack_func(&m, &n, x, fvec, iwa); }
    //
    // calc chisq
    for(int i=0; i<m; i++) { chisq[i] = fvec[i]*fvec[i]; }
    //
    // save to class
    this->chi2.resize(m);
    this->chi2.assign(chisq,chisq+m); //<TODO> need to check whether x and aCOE are identical.
    //
    // clean
    delete[] iwa;
    delete[] wa;
    delete[] fvec;
    // 
    // print on screen
    double sumchisq = 0.0;
    int debug = Input_debug;
    if(debug) {
        std::cout << "michi2MinPack:" << setw(15) << "fOBS" << setw(15) << "eOBS" << setw(15) << "fLIB1" << setw(15) << "fLIB2" << setw(15) << "fLIB3" << setw(15) << "fDIFF" << setw(15) << "chi^2" << std::endl;
        for(int i=0; i<m; i++) {
            // std::cout << "   " << setw(15) << fOBS[i] << setw(15) << eOBS[i] << setw(15) << fLIB[0][i] << std::endl ; // michi2
            // std::cout << "   " << setw(15) << fOBS[i] << setw(15) << eOBS[i] << setw(15) << fLIB[0][i] << setw(15) << fLIB[1][i] << std::endl ; // m2chi2
            std::cout << "              " << setw(15) << fOBS[i] << setw(15) << eOBS[i] << setw(15) << fLIB[0][i] << setw(15) << fLIB[1][i] << setw(15) << fLIB[2][i] << setw(15) << fOBS[i]-fLIB[0][i]*aCOE[0]-fLIB[1][i]*aCOE[1]-fLIB[2][i]*aCOE[2] << setw(15) << chi2[i] << "   + " << std::endl ; // m3chi2
        }
        std::cout << "michi2MinPack:" << " ncount=" << michi2MinPack_ncount << " info=" << info;
        std::cout << " n=" << n; for(int i=0; i<n; i++) { std::cout << " a[" << i << "]=" << aCOE[i]; }
        std::cout << " m=" << m; // for(int i=0; i<n; i++) { std::cout << " x[" << i << "]=" << x[i]; }
        for(int i=0; i<m; i++) { sumchisq+=chisq[i]; }
        std::cout << " sumchi2=" << sumchisq << std::endl;
    }
    //
    //
    switch (info) {
        case 0:
            std::cout << "michi2MinPack::fit() improper input parameters!" << std::endl;
            std::cout << "michi2MinPack:" << setw(15) << "fOBS" << setw(15) << "eOBS" << setw(15) << "fLIB1" << setw(15) << "fLIB2" << setw(15) << "fLIB3" << setw(15) << "fDIFF" << setw(15) << "chi^2" << std::endl;
            for(int i=0; i<m; i++) {
                std::cout << setw(14) << " " << setw(15) << fOBS[i] << setw(15) << eOBS[i] << setw(15) << fLIB[0][i] << setw(15) << fLIB[1][i] << setw(15) << fLIB[2][i] << setw(15) << fOBS[i]-fLIB[0][i]*aCOE[0]-fLIB[1][i]*aCOE[1]-fLIB[2][i]*aCOE[2] << setw(15) << chi2[i] << "   + " << std::endl ; // m3chi2
            }
            std::cout << "michi2MinPack:" << " ncount=" << michi2MinPack_ncount << " info=" << info;
            std::cout << " n=" << n; for(int i=0; i<n; i++) { std::cout << " a[" << i << "]=" << aCOE[i]; }
            std::cout << " m=" << m; // for(int i=0; i<n; i++) { std::cout << " x[" << i << "]=" << x[i]; }
            std::cout << " sumchi2=" << sumchisq << std::endl;
            break;
        case 1:
            // std::cout << "michi2MinPack::fit() algorithm estimates that the relative error in the sum of squares is at most tol." << std::endl;
            break;
        case 2:
            // std::cout << "michi2MinPack::fit() algorithm estimates that the relative error between x and the solution is at most tol." << std::endl;
            break;
        case 3:
            // std::cout << "michi2MinPack::fit() algorithm estimates that the relative error is at most tol." << std::endl;
            break;
        case 4:
            std::cout << "michi2MinPack::fit() fvec is orthogonal to the columns of the Jacobian to machine precision." << std::endl;
            break;
        case 5:
            std::cout << "michi2MinPack::fit() number of calls to fcn has reached or exceeded 200*(n+1)." << std::endl;
            break;
        case 6:
            std::cout << "michi2MinPack::fit() tol is too small. no further reduction in the sum of squares is possible." << std::endl;
            break;
        case 7:
            std::cout << "tol is too small. no further improvement in the approximate solution x is possible." << std::endl;
            break;
        default:
            break;
    }
}



double michi2MinPack::mean(std::vector<double> data)
{
    // return the mean value of *data
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / (double)data.size();
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / data.size() - mean * mean);
    return mean;
}





#endif
