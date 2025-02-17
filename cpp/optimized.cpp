#include <iostream>
#include <vector>
#include <numeric> 
#include <array>
#include <cmath>
#include <map>
#include <algorithm>
#include </home/zcandels/num_analysis_packages/eigen-3.4.0/Eigen/Dense>
#include "/home/zcandels/num_analysis_packages/nanoflann/include/nanoflann.hpp"
#include <unordered_map>
#include <valarray>
#include "/home/zcandels/geom/cpp/boundarynode.hpp"
#include <time.h>



double signed_distance(Eigen::Vector2d P, Eigen::Vector2d x0, Eigen::Vector2d u)
{
    Eigen::Vector2d v_perp(2);
    v_perp << -u[1], u[0];
    Eigen::Vector2d dist_vec = P - x0;

    //std::vector<double> P(std::begin(P), std::end())
    double dot_prod = dist_vec.dot(v_perp);

    return dot_prod;

}

std::vector<Eigen::Vector2d> find_closest_points(std::vector< Eigen::Vector2d > bdy_curve_pts, Eigen::Vector2d x0, Eigen::Vector2d unit_vec)
{
    double min_pos_dist = 100.;
    double min_neg_dist = 100.;

    Eigen::Vector2d closest_pos;
    Eigen::Vector2d closest_neg;
    std::vector<Eigen::Vector2d> closest_pts;

    for(int idx = 0; idx < bdy_curve_pts.size(); idx++)
    {
        Eigen::Vector2d point = bdy_curve_pts[idx];
        double dist = signed_distance(point, x0, unit_vec);

        if( dist > 0 && dist < min_pos_dist)
        {
            Eigen::Vector2d closest_pos = point;
            double min_pos_dist = dist;
        }
        else if (dist < 0 && std::abs(dist) < min_neg_dist)
        {
            Eigen::Vector2d closest_neg = point;
            double min_neg_dist = std::abs(dist);
        }

    }
    closest_pts.push_back(closest_neg);
    closest_pts.push_back(closest_pos);

    return closest_pts;
}


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


std::vector< Eigen::Vector2d > define_bdy_curve(int N_bdy_pts, double x_min, double x_max)
{
    double eps = 0.2;
    std::vector<double> x_vals = linspace(x_min - eps, x_max +eps, N_bdy_pts);
    std::vector<double> y_vals;

    std::vector< Eigen::Vector2d > bdy_curve_pts;

    for (int i = 0; i < N_bdy_pts; i++)
    {
        y_vals.push_back( pow(sin(5 * M_PI * x_vals[i]), 2.) );
        bdy_curve_pts.push_back( { x_vals[i], y_vals[i] } );

    }

    return bdy_curve_pts;
}


std::vector<std::vector<double>> label_solid_pts(std::vector<std::vector<double>>& grid_pts, 
                                                 const std::vector< Eigen::Vector2d >& bdy_curve_pts) 
{
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


std::vector<BoundaryNode> label_bdy_points(std::vector<std::vector<double>>& grid_pts,
    const std::vector<std::vector<double>>& solid_points) 
{
    std::vector<BoundaryNode> boundary_nodes;
    double eps = 1e-5;

    double dx = grid_pts[1][0] - grid_pts[0][0];
    double dy = grid_pts[1][1] - grid_pts[0][1];

    for (auto& grid_pt : grid_pts) 
    {
        if (grid_pt[2] == 2) continue; // Skip solid points

        Eigen::Vector2d P(grid_pt[0], grid_pt[1]);

        for (const auto& solid_pt : solid_points) 
        {
            Eigen::Vector2d Q(solid_pt[0], solid_pt[1]);
            Eigen::Vector2d Z = Q - P;

            if (Z.norm() <= sqrt(2.0) * std::max(dx, dy) + eps) 
            {
                grid_pt[2] = 1; // Mark as boundary point
                BoundaryNode x_b(P[0], P[1]);
                x_b.add_velocity_vec(Z.normalized());

                for (const auto& other_solid_pt : solid_points) 
                {
                    Eigen::Vector2d Q_other(other_solid_pt[0], other_solid_pt[1]);
                    if ((Q_other - P).norm() <= sqrt(2) * std::max(dx, dy) + eps) 
                    {
                        x_b.add_velocity_vec((Q_other - P).normalized());
                    }
                }

                boundary_nodes.push_back(x_b);
                break; // No need to check further
            }
        }
    }
    return boundary_nodes;
}


void distances_and_normals(std::vector<BoundaryNode>& boundary_nodes, std::vector< Eigen::Vector2d >& bdy_curve_pts)
{
    for(int j = 0; j < boundary_nodes.size(); j++)
    {
        Eigen::Vector2d x_b = boundary_nodes[j].getPosition();

        std::vector< Eigen::Vector2d > velocity_vecs = boundary_nodes[j].getVelDirections();

        for(int i = 0; i < velocity_vecs.size(); i++)
        {
            Eigen::Vector2d unit_vel_vec = boundary_nodes[j].getVelDirections()[i];
            std::vector<Eigen::Vector2d> closest_pts = find_closest_points(bdy_curve_pts, x_b, unit_vel_vec);
            Eigen::Vector2d pt1 = closest_pts[0];
            Eigen::Vector2d v = closest_pts[1] - closest_pts[0];
            v = v/v.norm();

            Eigen::Matrix2d intersection_mat;
            intersection_mat << v[0], -unit_vel_vec[0],
                                v[1], -unit_vel_vec[1];

            Eigen::Vector2d intersection_rhs_vec;
            intersection_rhs_vec << x_b[0] - pt1[0], x_b[1] - pt1[1];

            Eigen::Vector2d intersection_distances = intersection_mat.colPivHouseholderQr().solve(intersection_rhs_vec);

            double delta = intersection_distances[1];

            if( delta > 1)
                continue;
            else
            {
                boundary_nodes[j].setDistance(delta);
                Eigen::Vector2d normal_vec(-1, -v[0]/v[1]);
                normal_vec = normal_vec/normal_vec.norm();
                boundary_nodes[j].setNormals(normal_vec);
            }

        }

    }
}

void write_data(std::vector<std::vector<double>>& grid_pts, std::vector< Eigen::Vector2d >& bdy_curve_pts,
    std::vector<std::vector<double>>& solid_points, std::vector<BoundaryNode>& boundary_nodes )
{
    const int width = 10;
    const int precision = 5;

    FILE* outfile1 = fopen("boundary_pts.dat", "w");
    std::vector<std::vector<double>> boundary_node_positions;
    for(int i = 0; i < boundary_nodes.size(); i++)
    {
        Eigen::Vector2d x_b = boundary_nodes[i].getPosition();
        fprintf(outfile1, "%*.*f, %*.*f\n", width, precision, x_b[0], width, precision, x_b[1]);

    }
    fclose(outfile1);

    FILE* outfile2 = fopen("solid_pts.dat", "w");
    for(int i = 0; i < solid_points.size(); i++)
    {
        std::vector<double> x_s = solid_points[i];
        fprintf(outfile2, "%*.*f, %*.*f\n", width, precision, x_s[0], width, precision, x_s[1]);
    }
    fclose(outfile2);

    FILE* outfile3 = fopen("bdy_curve_pts.dat", "w");
    for(int i = 0; i < bdy_curve_pts.size(); i++)
    {
        Eigen::Vector2d x_c = bdy_curve_pts[i];
        fprintf(outfile2, "%*.*f, %*.*f\n", width, precision, x_c[0], width, precision, x_c[1]);
    }
    fclose(outfile3);

    FILE* outfile4 = fopen("grid_pts.dat", "w");
    for(int i = 0; i < grid_pts.size(); i++)
    {
        std::vector<double> x_h = grid_pts[i];
        fprintf(outfile2, "%*.*f, %*.*f\n", width, precision, x_h[0], width, precision, x_h[1]);
    }
    fclose(outfile4);
}



int main()
{

    std::clock_t start = clock(); 

    double alpha = 0.2;
    double R = 1;
    int N_bdy_pts = 3200;
    int N_repeats = 1;
    int n_grid_pts = 100;
    double x_min= -2;
    double x_max = 2;
    double y_min = 0;
    double y_max = 4;


    std::vector<std::vector<double>> grid_pts = make_grid(x_min, x_max, y_min, y_max, n_grid_pts);

    std::vector< Eigen::Vector2d > bdy_curve_pts = define_bdy_curve(N_bdy_pts, x_min, x_max);

    std::vector<std::vector<double>> solid_points = label_solid_pts(grid_pts, bdy_curve_pts);

    std::vector<BoundaryNode> boundary_nodes = label_bdy_points(grid_pts, solid_points);

    distances_and_normals(boundary_nodes, bdy_curve_pts);

    write_data(grid_pts, bdy_curve_pts, solid_points, boundary_nodes );

    std::cout << "Time taken: " << (double)(clock() - start)/CLOCKS_PER_SEC << std::endl;

}