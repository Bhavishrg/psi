#ifndef DISTANCE_H
#define DISTANCE_H

#include <Eigen/Dense>

double dist(const Eigen::VectorXd& a,
            const Eigen::VectorXd& b,
            int p = 2);

#endif