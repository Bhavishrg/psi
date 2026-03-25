#ifndef DATASET_H
#define DATASET_H

#include <vector>
#include <string>
#include <utility>
#include <Eigen/Dense>

using Vector = Eigen::VectorXd;

std::pair<std::vector<Vector>, std::vector<Vector>>
generate_sa_psi_dataset(
    int d,
    int N,
    int M,
    double delta,
    int U,
    std::string mode,
    int far_factor,
    int p
);

#endif