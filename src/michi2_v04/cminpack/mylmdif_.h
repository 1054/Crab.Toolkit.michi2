extern "C" {
    
typedef void (*cminpack_func_mn) (int m, int n, double *x, double *fvec, int iflag );

void lmdif (cminpack_func_mn fcn_mn, const int *m, const int *n, double *x, double *fvec, const double *ftol, const double *xtol, const double *gtol, const int *maxfev, const double *epsfcn, double *diag, const int *mode, const double *factor, const int *nprint, int *info, int *nfev, double *fjac, const int *ldfjac, int *ipvt, double *qtf, double *wa1, double *wa2, double *wa3, double *wa4);

}