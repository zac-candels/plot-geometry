#ifndef BOUNDARYNODEHEADERDEF
#define BOUNDARYNODEHEADERDEF

#include <vector>
#include <array>
#include <cmath>
#include <valarray>


class BoundaryNode
{
private:
    double x_position;
    double y_position;
    std::vector< std::valarray<double> > velocity_vecs;
    std::vector<double> distances;
    std::vector< std::array<double, 2> > normals;

public:
    BoundaryNode();
    BoundaryNode(double x, double y);

    void setDistance(double dist);
    void setNormals(std::array<double,2> normal);
    void add_velocity_vec(std::valarray<double> vel_vec);

    std::array<double, 2> getPosition();
    std::vector<double> getDistances();
    std::vector< std::array<double, 2> > getNormals();
    std::vector< std::valarray<double> > getVelDirections();
    

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

void BoundaryNode::setNormals(std::array<double, 2> normal)
{
    normals.push_back(normal);
}

void BoundaryNode::add_velocity_vec(std::valarray<double> vel_vec)
{
    velocity_vecs.push_back( vel_vec );
}

std::array<double, 2> BoundaryNode::getPosition()
{
    return {x_position, y_position};
}

std::vector< std::array<double, 2> > BoundaryNode::getNormals()
{
    return normals;
}

std::vector< std::valarray<double> > BoundaryNode::getVelDirections()
{
    return velocity_vecs;
}

#endif