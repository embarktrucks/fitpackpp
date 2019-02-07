#include <gtest/gtest.h>
#include <Eigen/Dense>
//#include "fitpackpp/BSplineCurve.h"
#include "fitpackpp/BSplineParCurve.h"
#include <iostream>

TEST(Fitpack,curfit) {
    //Create L-shaped data
    Eigen::MatrixXd pts = Eigen::MatrixXd::Zero(40,2);
    pts.col(0).head(20).setLinSpaced(20,0.0,10.0);
    pts.col(0).tail(20) = 10.0*Eigen::VectorXd::Ones(20);
    pts.col(1).tail(20).setLinSpaced(20,0.5,10.0);

    //fit spline to this, takes x,y as std::vector<double>
    std::vector<double> x(40);
    std::vector<double> y(40);
    for(size_t i=0; i<40; ++i){
        x.at(i) = pts(i,0);
        y.at(i) = pts(i,1);
    }
    Eigen::VectorXd s = Eigen::VectorXd::Zero(40);
    for(size_t i=1; i<40; ++i){
        s(i) = s(i-1)+ (pts.row(i)-pts.row(i-1)).norm();
    }
    std::vector<double> t(40);


    //this is only a function (x->y)!
    //fitpackpp::BSplineCurve bSpline(x,y);
    //Evaluate spline
    //bSpline.eval()

    //We want to have paramteric curve! (t->x,y)
    fitpackpp::BSplineParCurve bSpline(x,y,t);
    Eigen::VectorXd t_eval = Eigen::VectorXd::LinSpaced(5,0.3,0.8);
    Eigen::MatrixXd xy_interp = bSpline.eval(t_eval);
    std::cout << "xy_interp = \n" << xy_interp << "\n";

    Eigen::MatrixXd ders = bSpline.der(t_eval(0));
    std::cout << "ders = \n" << ders << "\n";

    ASSERT_TRUE(xy_interp.rows()==t_eval.rows());
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
