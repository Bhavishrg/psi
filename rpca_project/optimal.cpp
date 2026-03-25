#include "optimal.h"
#include <set>
#include <algorithm>

std::tuple<double,double,double>
optimal_b(
    const std::vector<Vector>& points,
    const Vector& v1,
    double delta,
    double epsilon
)
{
    std::vector<double> p;

    for(auto &x:points)
        p.push_back(x.dot(v1));

    std::sort(p.begin(),p.end());

    int N=p.size();

    double m=p[N/2];

    double left=m-epsilon;
    double right=m+epsilon;

    std::set<double> candidates={left,right,m};

    for(double pi:p)
    {
        double a=pi-delta;
        double b=pi+delta;

        if(left<=a && a<=right)
            candidates.insert(a);

        if(left<=b && b<=right)
            candidates.insert(b);
    }

    auto Z=[&](double b)
    {
        double s=0;

        for(double pi:p)
            if((pi-delta)<b && b<(pi+delta))
                s+=std::max(0.0,delta-abs(b-pi));

        return s;
    };

    double best_b=0;
    double best_val=1e18;

    for(double b:candidates)
    {
        double val=Z(b);

        if(val<best_val)
        {
            best_val=val;
            best_b=b;
        }
    }

    return {best_b,best_val,m};
}