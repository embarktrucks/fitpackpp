#ifndef BSPLINEPAR_CURVE_H
#define BSPLINEPAR_CURVE_H

#ifndef alloca
#define alloca __builtin_alloca
#endif

#include <vector>
#include <string>
#include <Eigen/Dense>

namespace fitpackpp
{

/**
 * @brief B-Spline curve
 * @details A wrapper for the curfit(), splev(), and splder() routines from FITPACK by  P. Dierckx:
 * http://www.netlib.org/fitpack/index.html
 * More specifically, we use the double-precision FITPACK version included with scipy:
 * https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 */
class BSplineParCurve
{
public:
    BSplineParCurve(std::vector<double> &x, std::vector<double> &y, std::vector<double> &t, int preferredDegree=3, double smoothing=0.0);

    ~BSplineParCurve();

    std::vector<double> knotX();
    std::vector<double> coefs();
    int degree();

    Eigen::MatrixXd eval(Eigen::VectorXd &t_eval);
    Eigen::MatrixXd der(Eigen::VectorXd &t_eval);
    Eigen::MatrixXd der(double t_eval);

private:
    int     k; // Spline degree
    int     n; // Number of knots
    double *t; // Knot coordinates
    double *c; // Spline coefficients
};

}

#endif
