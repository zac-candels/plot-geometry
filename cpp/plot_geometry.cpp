#include <iostream>
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

// End of the BoundaryNode class. Now we can move on to the other functions.



std::vector<std::vector<double>> make_grid(double x_min, double x_max, double y_min, double y_max, int n_grid_pts) {
    std::vector<double> x_pts(n_grid_pts);
    std::vector<double> y_pts(n_grid_pts);
    std::vector<std::vector<double>> grid_pts(n_grid_pts * n_grid_pts, std::vector<double>(3, 0.0));

    // Generate x and y coordinates
    for (int i = 0; i < n_grid_pts; ++i) {
        x_pts[i] = x_min + i * (x_max - x_min) / (n_grid_pts - 1);
        y_pts[i] = y_min + i * (y_max - y_min) / (n_grid_pts - 1);
    }

    // Assign coordinates to grid
    for (int i = 0; i < n_grid_pts; ++i) {
        for (int j = 0; j < n_grid_pts; ++j) {
            grid_pts[i + j * n_grid_pts][0] = x_pts[i];
            grid_pts[i + j * n_grid_pts][1] = y_pts[j];
        }
    }

    return grid_pts;
}


int main()
{
    double alpha = 0.2;
    double R = 1;
    int N_bdy_pts = 3200;
    int N_repeats = 1;
    int n_grid_pts = 10;
    double x_min= -2;
    double x_max = 2;
    double y_min = 0;
    double y_max = 4;

    double x = 1.8;
    double y = -3.5;
    std::vector< std::array<double, 2> > vel_dirs = { {2.3, -1.3}, {4.5, 0.1} };
    BoundaryNode x_b = BoundaryNode(2.3, 4.5, vel_dirs);

    std::array<double, 2> x_p = x_b.getPosition();

    std::vector<std::array<double, 2> > c_q = x_b.getVelDirections();



    std::cout << c_q[0][0] << ", " << c_q[0][1] << std::endl;


    std::vector<std::vector<double>> grid_pts = make_grid(x_min, x_max, y_min, y_max, n_grid_pts);



}