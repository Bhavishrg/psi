#include <iostream>

#include "dataset.h"
#include "rpca_tree.h"
#include "distance.h"

int main()
{
    int d=3;
    int k;
std::cout << "Enter k: ";
std::cin >> k;
std::vector<Node*> trees;


    auto data=generate_sa_psi_dataset(
        d,
        700,
        10000,
        5,
        255,
        "disjoint",
        3,
        2
    );

    auto centers=data.first;
    auto bob_points=data.second;

    double delta = 5;
    std::vector<Vector> Rql;   // bob points inside balls
std::vector<Vector> Rc;    // corresponding centers
int Rq = 0;

for(auto &b : bob_points)
{
    for(auto &a : centers)
    {
        if(dist(a,b,2) <= delta)
        {
            Rq++;
            Rql.push_back(b);
            Rc.push_back(a);
            break;
        }
    }
}

std::cout << "Actual matches: " << Rq << std::endl;


    for(int i = 0; i < k; i++)
{
    Matrix Rot = random_rotation_matrix(d);
    Node* T = R1pca_tree(centers,delta,Rot,0);
    trees.push_back(T);
}

int match_count = 0;

for(auto &q : bob_points)
{
    int r = 0;

    for(auto T : trees)
    {
        auto result = Rpca_search(T,q,delta,2);
        r = std::max(r,result.first);
    }

    if(r == 1)
        match_count++;
}

std::cout << "Protocol Matches: " << match_count << std::endl;
return 0;
}