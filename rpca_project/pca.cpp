#include "pca.h"

PCAResult my_pca(const std::vector<Vector>& X)
{
    int n=X.size();
    int d=X[0].size();

    Matrix M(n,d);

    for(int i=0;i<n;i++)
        M.row(i)=X[i];

    Vector mean=M.colwise().mean();
    M.rowwise()-=mean.transpose();

    Matrix Cov=(M.transpose()*M)/n;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(Cov);

    return {mean,Cov,solver.eigenvalues(),solver.eigenvectors()};
}