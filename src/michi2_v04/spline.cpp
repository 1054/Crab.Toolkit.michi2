#include <stdio.h>
#include <math.h>
#include <vector>
#include <iomanip>

using namespace std;

// https://stackoverflow.com/questions/1204553/are-there-any-good-libraries-for-solving-cubic-splines-in-c

// https://en.wikipedia.org/w/index.php?title=Spline_%28mathematics%29&oldid=288288033#Algorithm_for_computing_natural_cubic_splines


/*
    A common spline is the natural cubic spline of degree 3 with continuity C2. The word "natural" means that the second derivatives of the spline polynomials are set equal to zero at the endpoints of the interval of interpolation

    {\displaystyle S''(a)\,=S''(b)=0.} {\displaystyle S''(a)\,=S''(b)=0.}
    This forces the spline to be a straight line outside of the interval, while not disrupting its smoothness.
 
    Algorithm for computing natural cubic splines
    Cubic splines are of the form
    S_j (x) = a_j + b_j * (x - x_j) + c_j * (x - x_j)^2 + d_j * (x - x_j)^3
    Given set of coordinates (x_i, y_i), i=1,n, we wish to find set of n splines S_i (x) for i=1,n

 */



typedef struct {
    double a;
    double b;
    double c;
    double d;
    double x;
} SplineSet;



vector<double> spline(vector<double> &x, vector<double> &y, vector<double> &output_x)
{
    // x must be sorted
    
    
    int i=0, j=0, k=0;
    int n = x.size()-1;
    vector<double> a(n+1);
    vector<double> b(n);
    vector<double> d(n);
    vector<double> h;
    vector<double> alpha;
    
    // set a_i = y_i for i = 0,n
    a.insert(a.begin(), y.begin(), y.end());
    
    // set h_i = (x_i+1 - x_i) for i = 0,n-1
    for(i=0; i<n; i++)
        h.push_back(x[i+1]-x[i]);
    
    // set alpha_i = 3/h_i * (a_i+1 - a_i) - ... for i=1,n-1
    for(i=1; i<n; i++)
        alpha.push_back(3.0/h[i]*(a[i+1]-a[i])-3.0/h[i-1]*(a[i]-a[i-1]));
    
    vector<double> c(n+1); c[n]=0;
    vector<double> l(n+1); l[0]=1; l[n]=1;
    vector<double> mu(n+1); mu[0]=0;
    vector<double> z(n+1); z[0]=0; z[n]=0; // the second derivative at the ith knot
    
    for(i=1; i<n; i++) {
        l[i] = 2.0*(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }
    
    for(j=n-1; j>=0; j--) {
        c[j] = z[j] - mu[j] * c[j+1];
        b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3.0;
        d[j] = (c[j+1]-c[j])/(3.0*h[j]);
    }
    
    vector<SplineSet> output_set(n);
    for(i=0; i<n; i++) {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = c[i];
        output_set[i].d = d[i];
        output_set[i].x = x[i];
    }

    
    //vector<double> output_y(n+1); // SUPER BUG!!! WAISTED ME 5 HOURS! FROM 2017-10-01 23h59m TO 2017-10-02 05h12m!!! WTF!!! 
    vector<double> output_y(output_x.size());
    
    for(k=0; k<output_x.size(); k++) {
        output_y[k] = 0;
        if(output_x.at(k)>=x.front() && output_x.at(k)<=x.back()) {
            double min_val = -99;
            int min_ind = -1;
            for(i=0; i<n; i++) {
                if(min_val<0 || min_val>fabs(output_x.at(k)-output_set[i].x)) {
                    min_val = fabs(output_x.at(k)-output_set[i].x);
                    min_ind = i;
                }
            }
            i=min_ind;
            output_y[k] = output_set[i].a +
                          output_set[i].b * (output_x.at(k) - output_set[i].x) +
                          output_set[i].c * (output_x.at(k) - output_set[i].x) * (output_x.at(k) - output_set[i].x) +
                          output_set[i].d * (output_x.at(k) - output_set[i].x) * (output_x.at(k) - output_set[i].x) * (output_x.at(k) - output_set[i].x);
        }
    }
    
//    std::cout << "spline: debugging: x=0x" << std::hex << (size_t)&x << std::endl;
//    std::cout << "spline: debugging: y=0x" << std::hex << (size_t)&y << std::endl;
//    std::cout << "spline: debugging: output_x=0x" << std::hex << (size_t)&output_x << std::endl;
//    std::cout << "spline: debugging: a=0x" << std::hex << (size_t)&a << std::endl;
//    std::cout << "spline: debugging: b=0x" << std::hex << (size_t)&b << std::endl;
//    std::cout << "spline: debugging: d=0x" << std::hex << (size_t)&d << std::endl;
//    std::cout << "spline: debugging: h=0x" << std::hex << (size_t)&h << std::endl;
//    std::cout << "spline: debugging: output_set=0x" << std::hex << (size_t)&output_set << std::endl;
//    std::cout << "spline: debugging: output_y=0x" << std::hex << (size_t)&output_y << std::endl;
    
    return output_y;
}





