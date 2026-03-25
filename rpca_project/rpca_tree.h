#ifndef RPCA_TREE_H
#define RPCA_TREE_H

#include <Eigen/Dense>
#include <vector>

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

struct Node
{
    bool leaf;

    Vector point;

    Vector v;
    double threshold;

    Node* left;
    Node* right;

    Node(Vector p);
    Node(Vector v,double t,Node* l,Node* r);
};

Matrix random_rotation_matrix(int d);

Node* R1pca_tree(
    std::vector<Vector> points,
    double delta,
    Matrix Rot,
    int i
);

std::pair<int,Vector>
Rpca_search(
    Node* T,
    const Vector& q,
    double delta,
    int p
);

#endif