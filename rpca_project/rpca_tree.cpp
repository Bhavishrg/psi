#include "rpca_tree.h"
#include "distance.h"
#include "pca.h"
#include "optimal.h"

#include <random>
#include <vector>
#include <Eigen/Dense>

Node::Node(Vector p)
{
    leaf=true;
    point=p;
    left=nullptr;
    right=nullptr;
}

Node::Node(Vector v_,double t,Node* l,Node* r)
{
    leaf=false;
    v=v_;
    threshold=t;
    left=l;
    right=r;
}

Matrix random_rotation_matrix(int d)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> g(0,1);

    Matrix A(d,d);

    for(int i=0;i<d;i++)
        for(int j=0;j<d;j++)
            A(i,j)=g(gen);

    Eigen::HouseholderQR<Matrix> qr(A);

    Matrix Q=qr.householderQ();

    if(Q.determinant()<0)
        Q.col(0)=-Q.col(0);

    return Q;
}

Node* R1pca_tree(
    std::vector<Vector> points,
    double delta,
    Matrix Rot,
    int i
)
{
    if(points.size()<=1)
        return new Node(points[0]);

    auto pca=my_pca(points);

    Vector v1=Rot*pca.eigenvectors.col(i);

    auto [tc,val,m]=optimal_b(points,v1,delta,delta);

    std::vector<Vector> left,right;

    for(auto &p:points)
    {
        if(p.dot(v1)<tc)
            left.push_back(p);
        else
            right.push_back(p);
    }

    if(left.empty()||right.empty())
    {
        int mid=points.size()/2;

        left.assign(points.begin(),points.begin()+mid);
        right.assign(points.begin()+mid,points.end());
    }

    Node* L=R1pca_tree(left,delta,Rot,i);
    Node* R=R1pca_tree(right,delta,Rot,i);

    return new Node(v1,tc,L,R);
}

std::pair<int,Vector>
Rpca_search(
    Node* T,
    const Vector& q,
    double delta,
    int p
)
{
    while(!T->leaf)
    {
        if(T->v.dot(q)<=T->threshold)
            T=T->left;
        else
            T=T->right;
    }

    Vector x=T->point;

    int result=dist(q,x,p)<=delta;

    return {result,x};
}