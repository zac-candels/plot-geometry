#ifndef BOUNDARYNODEHEADERDEF
#define BOUNDARYNODEHEADERDEF

#include <vector>
#include <array>
#include <cmath>


class BoundaryNode
{
private:
    double x_position;
    double y_position;
    std::vector< std::array<double,2> > velocity_vecs;
    std::vector<double> distances;
    std::vector< std::array<double, 2> > normals;

public:
    BoundaryNode();
    BoundaryNode(double x, double y, std::vector< std::array<double,2> > directions);

    void setDistance(double dist);
    void setNormals(std::array<double,2> normal);

    std::array<double, 2> getPosition();
    std::vector<double> getDistances();
    std::vector< std::array<double, 2> > getNormals();
    std::vector< std::array<double, 2> > getVelDirections();
    

};


BoundaryNode::BoundaryNode(double x, double y, std::vector< std::array<double,2> > directions)
{
    x_position = x;
    y_position = y;
    velocity_vecs = directions;

}

void BoundaryNode::setDistance(double dist)
{
    distances.push_back(dist);
}

void BoundaryNode::setNormals(std::array<double, 2> normal)
{
    normals.push_back(normal);
}

std::array<double, 2> BoundaryNode::getPosition()
{
    return {x_position, y_position};
}

std::vector< std::array<double, 2> > BoundaryNode::getNormals()
{
    return normals;
}

std::vector< std::array<double, 2> > BoundaryNode::getVelDirections()
{
    return velocity_vecs;
}

#endif