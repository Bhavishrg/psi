#include "dataset.h"
#include "distance.h"
#include <random>

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
)
{
    std::mt19937 gen(42);
    std::uniform_int_distribution<> distU(0, U-1);

    std::vector<Vector> centers;
    double sep;

    if(mode=="disjoint") sep = 2*delta+1;
    else if(mode=="disjoint_and_far") sep = far_factor*delta;
    else sep = delta;

    while(centers.size()<N)
    {
        Vector c(d);
        for(int i=0;i<d;i++) c[i]=distU(gen);

        bool ok=true;
        for(auto &cp:centers)
            if(dist(c,cp,p)<sep) ok=false;

        if(ok) centers.push_back(c);
    }

    std::vector<Vector> bob_points;

    for(int i=0;i<M;i++)
    {
        Vector b(d);
        for(int j=0;j<d;j++) b[j]=distU(gen);
        bob_points.push_back(b);
    }

    return {centers,bob_points};
}