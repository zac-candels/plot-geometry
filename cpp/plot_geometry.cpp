#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include </home/zcandels/num_analysis_packages/eigen-3.4.0/Eigen/Dense>
#include <unordered_map>
#include "boundarynode.hpp"

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

// End of the BoundaryNode class. Now we can move on to the other functions.



std::vector<std::vector<double>> make_grid(double x_min, double x_max, double y_min, double y_max, int n_grid_pts) 
{
    std::vector<double> x_pts(n_grid_pts);
    std::vector<double> y_pts(n_grid_pts);
    std::vector<std::vector<double>> grid_pts(n_grid_pts * n_grid_pts, std::vector<double>(3, 0.0));

    // Generate x and y coordinates
    for (int i = 0; i < n_grid_pts; ++i) 
    {
        x_pts[i] = x_min + i * (x_max - x_min) / (n_grid_pts - 1);
        y_pts[i] = y_min + i * (y_max - y_min) / (n_grid_pts - 1);
    }

    // Assign coordinates to grid
    for (int i = 0; i < n_grid_pts; ++i) 
    {
        for (int j = 0; j < n_grid_pts; ++j) 
        {
            grid_pts[i + j * n_grid_pts][0] = x_pts[i];
            grid_pts[i + j * n_grid_pts][1] = y_pts[j];
        }
    }

    return grid_pts;
}


std::vector< std::array<double, 2> > define_bdy_curve(int N_bdy_pts, double x_min, double x_max)
{
    double eps = 0.2;
    std::vector<double> x_vals = linspace(x_min, x_max, N_bdy_pts);
    std::vector<double> y_vals;

    std::vector< std::array<double, 2> > bdy_curve_pts;

    for (int i = 0; i < N_bdy_pts; i++)
    {
        y_vals.push_back( pow(sin(5 * M_PI * x_vals[i]), 2) );
        bdy_curve_pts.push_back( { x_vals[i], y_vals[i] } );

    }

    return bdy_curve_pts;
}


std::vector<std::vector<double>> label_solid_pts(std::vector<std::vector<double>>& grid_pts, 
                                                 const std::vector< std::array<double, 2> >& bdy_curve_pts) {
    std::vector<std::vector<double>> solid_points;

    // Linear interpolation using Eigen
    Eigen::VectorXd x_vals(bdy_curve_pts.size()), y_vals(bdy_curve_pts.size());

    for (size_t i = 0; i < bdy_curve_pts.size(); ++i) 
    {
        x_vals[i] = bdy_curve_pts[i][0];
        y_vals[i] = bdy_curve_pts[i][1];
    }

    for (size_t i = 0; i < grid_pts.size(); ++i) 
    {
        double x_mesh_pt = grid_pts[i][0];
        double y_mesh_pt = grid_pts[i][1];

        // Find closest x values for interpolation
        auto it = std::lower_bound(x_vals.data(), x_vals.data() + x_vals.size(), x_mesh_pt);
        int idx = std::distance(x_vals.data(), it);
        if (idx == 0 || idx == x_vals.size()) continue; // Out of bounds

        // Linear interpolation
        double x1 = x_vals[idx - 1], x2 = x_vals[idx];
        double y1 = y_vals[idx - 1], y2 = y_vals[idx];
        double y_curve = y1 + (y2 - y1) * (x_mesh_pt - x1) / (x2 - x1);

        if (y_mesh_pt < y_curve) {
            grid_pts[i][2] = 2;
            solid_points.push_back({x_mesh_pt, y_mesh_pt});
        }
    }

    return solid_points;
}


/*
// BoundaryNode class definition
class BoundaryNode {
public:
    double x, y;
    std::vector<Eigen::Vector2d> velocity_vecs;
    std::vector<double> distances;
    std::vector<Eigen::Vector2d> normals;
    
    BoundaryNode(double x_, double y_, const std::vector<Eigen::Vector2d>& velocity_vecs_)
        : x(x_), y(y_), velocity_vecs(velocity_vecs_) {}
};
*/




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



    //std::cout << c_q[0][0] << ", " << c_q[0][1] << std::endl;


    std::vector<std::vector<double>> grid_pts = make_grid(x_min, x_max, y_min, y_max, n_grid_pts);

    std::vector< std::array<double, 2> > bdy_curve_pts = define_bdy_curve(N_bdy_pts, x_min, x_max);


    std::vector<std::vector<double>> solid_points = label_solid_pts(grid_pts, bdy_curve_pts);



}