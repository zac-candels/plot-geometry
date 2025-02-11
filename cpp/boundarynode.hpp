#ifndef BOUNDARYNODEHEADERDEF
#define BOUNDARYNODEHEADERDEF

#include <vector>
#include <array>
#include <cmath>
#include <valarray>
#include </home/zcandels/num_analysis_packages/eigen-3.4.0/Eigen/Dense>


class BoundaryNode
{
private:
    double x_position;
    double y_position;
    std::vector< Eigen::Vector2d > velocity_vecs;
    std::vector<double> distances;
    std::vector< Eigen::Vector2d > normals;

public:
    BoundaryNode();
    BoundaryNode(double x, double y);

    void setDistance(double dist);
    void setNormals(Eigen::Vector2d normal);
    void add_velocity_vec(Eigen::Vector2d vel_vec);

    Eigen::Vector2d getPosition();
    std::vector<double> getDistances();
    std::vector< Eigen::Vector2d > getNormals();
    std::vector< Eigen::Vector2d > getVelDirections();
    

};


BoundaryNode::BoundaryNode(double x, double y)
{
    x_position = x;
    y_position = y;

}

void BoundaryNode::setDistance(double dist)
{
    distances.push_back(dist);
}

void BoundaryNode::setNormals(Eigen::Vector2d normal)
{
    normals.push_back(normal);
}

void BoundaryNode::add_velocity_vec(Eigen::Vector2d vel_vec)
{
    velocity_vecs.push_back( vel_vec );
}

Eigen::Vector2d BoundaryNode::getPosition()
{
    return {x_position, y_position};
}

std::vector< Eigen::Vector2d > BoundaryNode::getNormals()
{
    return normals;
}

std::vector< Eigen::Vector2d > BoundaryNode::getVelDirections()
{
    return velocity_vecs;
}

#endif