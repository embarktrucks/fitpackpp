#include <assert.h>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <ios>

#include "BSplineParCurve.h"

#include "FCMangle.h"

using namespace fitpackpp;

extern "C" {
    //http://www.netlib.org/dierckx/parcur.f
    void parcur(int *iopt, int *ipar, int *idim, int *m,
                double *u, int *mx, double *x, double *w,
                double *ub, double *ue, int* k, double *s,
                int* nest, int* n, double * t, int * nc,
                double *c, double * fp,
                double* wrk, int* lwrk, int* iwrk, int* ier);

    void curev(int *idim, double *t, int* n, double *c, int *nc,
               int *k, double *u, int *m, double *x, int * mx, int * ier);

    void cualde(int *idim, double *t, int* n, double *c,int *nc,
                int *k1,double *u,double *d, int* nd, int* ier);


}

///TODO allow for more than 2d curve
BSplineParCurve::BSplineParCurve(std::vector<double> &x_data,
                                 std::vector<double> &y_data,
                                 std::vector<double> &t_data,
                                 int preferredDegree,
                                 double smoothing) {
    //in the order of parcur.f
    int iopt = 0;  // Compute a smoothing spline
    int ipar = 0;
    int idim = 2;
    int m = (int) x_data.size(); // Number of data points
    double * u = t_data.data();
    int mx = idim*m; //size of x array
    double * x = new double[m*idim];
    //for(size_t ii=0; ii<idim; ++ii){
    for(size_t ii=0; ii<m; ++ii){
        //x(idim*(i-1)+j) must contain the j-th coorinate of the i-th data
        //point for i=1,2,...,m and j=1,2,...,idim
        x[idim*ii+0] = x_data.at(ii);
        x[idim*ii+1] = y_data.at(ii);
    }
    //}
    double *w = new double[m];
    std::fill(w, w + m, 1.0);
    double ub=t_data.front();
    double ue=t_data.back();
    // The actual degree of the spline must be less than m
    k = preferredDegree;
    if (k >= m) {
        k = m - 1;
        std::cerr << "WARNING:  Too few data points (" << m << ") to create B-Spline curve of order " << preferredDegree << ". Reducing order to " << k << "." << std::endl;
    }
    double s = smoothing; //Smooting factor
    int nest = m + k + 1;               // Over-estimate the number of knots
    //total number of knots n
    t = new double[nest];               // Knots
    std::fill(t, t + nest, 0.0);
    int nc = nest*idim;
    c = new double[nest*idim];               // Coefficients
    double fp = 0.0; // Weighted sum of squared residuals

    // Allocate working memory required by curfit
    int     lwrk = m*(k+1)+nest*(6+idim+3*k);
    double *wrk  = new double[lwrk];
    std::fill(wrk, wrk + lwrk, 0.0);

    int    *iwrk = new int   [nest];
    std::fill(iwrk, iwrk + nest, 0);
    int ier = 0;

    parcur(&iopt,&ipar,&idim,&m,u,&mx,x,w,&ub,&ue,&k,&s,&nest,&n,t,&nc,c,&fp,wrk,&lwrk,iwrk,&ier);
    if (ier > 0) {
        if (ier >= 10) {
            std::stringstream s;
            s << "Error fitting B-Spline curve using parcur(): " << ier;
            throw std::runtime_error(s.str());
        } else {
            std::cerr << "WARNING:  Non-fatal error while fitting B-Spline curve using parcur(): " << ier << std::endl;
        }
    }

    // De-allocate temporary memory
    delete[] x;
    delete[] w;
    delete[] wrk;
    delete[] iwrk;
}

BSplineParCurve::~BSplineParCurve() {
    delete[] t;
    delete[] c;
}

Eigen::MatrixXd BSplineParCurve::eval(Eigen::VectorXd & t_eval) {
    //void splev (double *t, int *n, double *c, int *k,
    //            double *x, double *y, int *m, int *e, int *ier);

    int idim=2;
    //we have t from fitting (Knots)
    //and n (number of knots)
    //and c
    // coefficients of the spline sj(u) will be given
    // in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
    int nc = n-k-1;
    //We get errors with this, maybe the memory layout is bad for fortran!
    //double * u = t_eval.data();
    int m = t_eval.rows();
    double * u = new double[m];
    for(int ii=0; ii<m; ++ii) {
        u[ii] = t_eval(ii);
    }

    int mx = m*idim;
    double * x = new double[m*idim];

    int ier;

    //void curev(int *idim, double *t, int* n, double *c, int *nc,
    //               int *k, double *u, int *m, int *mx, double * x, int * ier);

    //t(k+1) <= u(i) <= u(i+1) <= t(n-k) , i=1,2,...,m-1.


    curev(&idim, t, &n, c, &nc,&k,u,&m,x,&mx,&ier);
    if (ier > 0) {
        if (ier >= 10) {
            std::stringstream s;
            s << "Error eval B-Spline curve using curev(): " << ier;
            throw std::runtime_error(s.str());
        } else {
            std::cerr << "WARNING:  Non-fatal error while eval B-Spline curve using curev(): " << ier << std::endl;
        }
    }


    Eigen::MatrixXd x_eval(m,idim);
    for(int i=0; i<m; ++i)
        for(int j=0; j<idim; ++j)
            x_eval(i,j) = x[idim*i+j];

    delete[] u;
    delete[] x;

    return x_eval;
}

Eigen::MatrixXd BSplineParCurve::der(double t_eval) {
    int idim=2;
    //we have t
    //have n
    //have c
    int nc = n-k-1;
    int k1 = k+1;
    int nd = k1*idim;
    double * d = new double[idim*k1];
    int ier;
    double u = 0.0;

    //Just first and second derivatives??


    u = t_eval;

    cualde(&idim,t,&n,c,&nc,&k1,&u,d,&nd,&ier);

    if (ier > 0) {
        if (ier >= 10) {
            std::stringstream s;
            s << "Error eval B-Spline curve using cualde(): " << ier;
            throw std::runtime_error(s.str());
        } else {
            std::cerr << "WARNING:  Non-fatal error while eval B-Spline curve using cualde(): " << ier << std::endl;
        }
    }

    //array,length nd,giving the different curve derivatives.
    // d(idim*l+j) will contain the j-th coordinate of the l-th
    // derivative of the curve at the point u.

    Eigen::MatrixXd ders(k1,idim);

    for(int l=0; l<k1; ++l){
        //std::cout << "d = ";
        for(int j=0; j<idim; ++j) {
            //std::cout << d[idim*l+j] << ",";
            ders(l,j) = d[idim*l+j];
        }
        //std::cout << "\n";
    }


    delete[] d;
    return ders;
}
Eigen::MatrixXd BSplineParCurve::der(Eigen::VectorXd &t_eval) {
    int idim=2;
    //we have t
    //have n
    //have c
    int nc = n-k-1;
    int k1 = k+1;
    int nd = k1*idim;
    double * d = new double[idim*k1];
    int ier;
    double u = 0.0;

    //Just first and second derivatives??

    for(int i=0; i<t_eval.rows(); ++i) {
        u = t_eval(i);

        cualde(&idim,t,&n,c,&nc,&k1,&u,d,&nd,&ier);

        if (ier > 0) {
            if (ier >= 10) {
                std::stringstream s;
                s << "Error eval B-Spline curve using cualde(): " << ier;
                throw std::runtime_error(s.str());
            } else {
                std::cerr << "WARNING:  Non-fatal error while eval B-Spline curve using cualde(): " << ier << std::endl;
            }
        }

        //array,length nd,giving the different curve derivatives.
        // d(idim*l+j) will contain the j-th coordinate of the l-th
        // derivative of the curve at the point u.
        for(int l=0; l<k1; ++l){
            std::cout << "d = ";
            for(int j=0; j<idim; ++j) {
                std::cout << d[idim*l+j] << ",";
            }
            std::cout << "\n";
        }
    }

    delete[] d;
    return Eigen::MatrixXd::Zero(5,5);
}
