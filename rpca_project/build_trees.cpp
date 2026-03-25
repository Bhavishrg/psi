#include <iostream>
#include <vector>

#include "rpca_tree.h"

int main()
{
    int d;
    int N;
    int k;

    std::cin >> d >> N >> k;

    std::vector<Vector> centers;

    for(int i=0;i<N;i++)
    {
        Vector v(d);

        for(int j=0;j<d;j++)
            std::cin >> v[j];

        centers.push_back(v);
    }

    double delta = 5;

    std::vector<Node*> trees;

    for(int i=0;i<k;i++)
    {
        Matrix Rot = random_rotation_matrix(d);
        Node* T = R1pca_tree(centers,delta,Rot,0);
        trees.push_back(T);

        std::cout << "Tree " << i << " built" << std::endl;
    }

    return 0;
}