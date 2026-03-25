#ifndef PCA_H
#define PCA_H

#include <Eigen/Dense>
#include <vector>

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

struct PCAResult
{
    Vector mean;
    Matrix covariance;
    Vector eigenvalues;
    Matrix eigenvectors;
};

PCAResult my_pca(const std::vector<Vector>& X);

#endif