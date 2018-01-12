/* 
    These codes are from "Numerical Recipes in C (Press et al. 1992).pdf"
    The Aim is to compute the $\Delta \chi^2$ for N free parameters of interest in least chi-square fitting.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip> // std::setprecision
#include "gamma.cpp"
#define ITMAX 100 /* Maximum allowed number of iterations. Relative accuracy. */
#define EPS 3.0e-7 /* Relative accuracy. */
#define FPMIN 1.0e-30 /* Number near the smallest representable doubleing-point number. */

using namespace std;


// nrerror()
// Numerical Recipes standard error handler
void nrerror(std::string error_text)
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text.c_str());
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}


// gcf()
// Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation as gammcf.
// Also returns lnΓ(a) as gln.
void gcf(double *gammcf, double a, double x, double *gln)
{
    double gammln(double xx);
    void nrerror(std::string error_text);
    int i;
    double an,b,c,d,del,h;
    *gln=gammln(a);
    b=x+1.0-a; // Set up for evaluating continued fraction by modified Lentz’s method (§5.2) with b0 = 0.
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=ITMAX;i++) { // Iterate to convergence.
        an = -i*(i-a);
        b += 2.0;
        d=an*d+b;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=b+an/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
    *gammcf=exp(-x+a*log(x)-(*gln))*h; // Put factors in front.
}


// gser()
// Returns the incomplete gamma function P (a, x) evaluated by its series representation as gamser.
// Also returns lnΓ(a) as gln.
void gser(double *gamser, double a, double x, double *gln)
{
    double gammln(double xx);
    void nrerror(std::string error_text);
    int n;
    double sum,del,ap;
    *gln=gammln(a); if (x <= 0.0) {
        if (x < 0.0) nrerror("x less than 0 in routine gser"); *gamser=0.0;
        return;
    } else { ap=a;
        del=sum=1.0/a;
        for (n=1;n<=ITMAX;n++) {
            ++ap;
            del *= x/ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS) {
                *gamser=sum*exp(-x+a*log(x)-(*gln));
                return; }
        }
        nrerror("a too large, ITMAX too small in routine gser"); return;
    }
}


// Returns the incomplete gamma function P(a, x).
double gammp(double a, double x)
{
    void gcf(double *gammcf, double a, double x, double *gln);
    void gser(double *gamser, double a, double x, double *gln);
    void nrerror(std::string error_text);
    double gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
    if (x < (a+1.0)) { // Use the series representation.
        gser(&gamser,a,x,&gln); return gamser;
    } else { // Use the continued fraction representation and take its complement.
        gcf(&gammcf,a,x,&gln); return 1.0-gammcf;
    }
}


// Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
double gammq(double a, double x)
{
    void gcf(double *gammcf, double a, double x, double *gln);
    void gser(double *gamser, double a, double x, double *gln);
    void nrerror(std::string error_text);
    double gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
    if (x < (a+1.0)) { // Use the series representation and take its complement.
        gser(&gamser,a,x,&gln); return 1.0-gamser;
    } else { // Use the continued fraction representation.
        gcf(&gammcf,a,x,&gln); return gammcf;
    }
}


// main()
// use the routine gammq and a simple root-finding routine (e.g., bisection) to find ∆ such that gammq(ν/2, ∆/2) = 1 − p.
int main(int argc, char **argv) {
    int n_freepar = 0; // 1; // 6;
    double Delta_chisq = 0.00; // 1.00; // 7.04;
    double p_confident = 0.68;
    double p_gamma_func = 0.00; // gammq(double(n_freepar)/2.0, Delta_chisq/2.0);
    // read user input
    //std::cout << argc << std::endl; // argv[0] is the program itself. User input starts from argv[1].
    if(argc>1) {
        n_freepar = std::stoi(argv[1]);
    } else {
        std::string str_n_freepar;
        std::cout << "Please input number of free parameter of interest: ";
        std::cin >> str_n_freepar;
        //std::cout << "User input: " << str_n_freepar << std::endl;
        n_freepar = std::stoi(str_n_freepar);
    }
    // check n_freepar
    if(n_freepar<=0) nrerror("Invalid number of free parameter! Must be a positive integer number!");
    std::cout << "ν = " << n_freepar << std::endl;
    // read user input of the desired p_confident
    if(argc>2) {
        p_confident = std::stod(argv[2]);
    }
    // check n_freepar
    if(p_confident<=0.0) nrerror("Invalid value of desired confidence! Must be a positive float number!");
    std::cout << "p = " << std::fixed << std::setprecision(6) << p_confident << std::endl;
    // now use bisection root-finding method to find the 'Delta_chisq' which makes '1-p_gamma_func == p_confident'.
    double t_bisec = 1.0e+30;
    double t_lower = 0.0; // possible lower value for Delta_chisq
    double t_upper = 1.0e+10; // possible upper value for Delta_chisq
    double t_middle = 0.0;
    double t_lower_gamma_func = gammq(double(n_freepar)/2.0, t_lower/2.0);
    double t_middle_gamma_func = 0.0;
    long t_count = 0;
    while (fabs(t_bisec) > FPMIN) {
        t_count++;
        t_middle = (t_lower+t_upper)/2.0;
        t_middle_gamma_func = gammq(double(n_freepar)/2.0, t_middle/2.0);
        t_bisec = (1-t_middle_gamma_func) - p_confident;
        //std::cout << std::scientific << std::setprecision(15) << t_middle << " " << t_upper-t_lower << " " << t_bisec << std::endl; // dump & debug
        if(fabs(t_bisec) > FPMIN) {
            //std::cout << std::scientific << std::setprecision(15) << (1-t_middle_gamma_func-p_confident)*(1-t_lower_gamma_func-p_confident) << std::endl;
            if((1-t_middle_gamma_func-p_confident)*(1-t_lower_gamma_func-p_confident)<0) {
                t_upper = t_middle;
            } else {
                t_lower = t_middle;
                t_lower_gamma_func = t_middle_gamma_func;
            }
        }
        if(fabs(t_upper-t_lower) < EPS) {
            break;
        }
        if(t_count>=100000) {
            std::cout << "Warning! Too many iterations! The result will likely be wrong!" << std::endl;
            break;
        }
    }
    // now found the 'Delta_chisq' which makes '1-p_gamma_func == p_confident'.
    Delta_chisq = t_middle;
    p_gamma_func = t_middle_gamma_func;
    // print
    //std::cout << "ν = " << std::fixed << std::setprecision(0) << n_freepar << std::endl;
    std::cout << "∆ = " << std::fixed << std::setprecision(6) << Delta_chisq << std::endl;
    //std::cout << "1 - gammq(ν/2, ∆/2) = " << std::fixed << std::setprecision(6) << 1 - p_gamma_func << std::endl;
    //std::cout << "p_confident = " << std::fixed << std::setprecision(6) << p_confident << std::endl;
    return 0;
}










