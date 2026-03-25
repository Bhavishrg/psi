#ifndef OPTIMAL_H
#define OPTIMAL_H

#include <Eigen/Dense>
#include <vector>
#include <tuple>

using Vector = Eigen::VectorXd;

std::tuple<double,double,double>
optimal_b(
    const std::vector<Vector>& points,
    const Vector& v1,
    double delta,
    double epsilon
);

#endif