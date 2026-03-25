#include "distance.h"
#include <stdexcept>
#include <Eigen/Dense>

double dist(const Eigen::VectorXd& a,
            const Eigen::VectorXd& b,
            int p)
{
    if (p == 1)
        return (a - b).lpNorm<1>();

    else if (p == 2)
        return (a - b).norm();

    else if (p == -1)
        return (a - b).lpNorm<Eigen::Infinity>();

    else
        throw std::invalid_argument("p must be 1,2,-1");
}