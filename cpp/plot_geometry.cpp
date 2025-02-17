#include <iostream>
#include <vector>
#include <numeric> 
#include <array>
#include <cmath>
#include <map>
#include <algorithm>
#include </home/zcandels/num_analysis_packages/eigen-3.4.0/Eigen/Dense>
#include <unordered_map>
#include <valarray>
#include <time.h>
#include <omp.h>

struct BdyNode
{
    double x_position;
    double y_position;
    Eigen::Vector2d Position;
    std::vector< Eigen::Vector2d > velocity_vecs;
    std::vector< std::vector<double> > vecs;
    std::vector<double> distances;
    std::vector< Eigen::Vector2d > normals;
    Eigen::MatrixXd matrix;
};


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

    #pragma omp parallel for
    for(int idx = 0; idx < bdy_curve_pts.size(); idx++)
    {
        Eigen::Vector2d point = bdy_curve_pts[idx];
        double dist = signed_distance(point, x0, unit_vec);


        #pragma omp critical
        {
            if( dist > 0 && dist < min_pos_dist)
            {
                closest_pos = point;
                min_pos_dist = dist;
            }
            else if (dist < 0 && std::abs(dist) < min_neg_dist)
            {
                closest_neg = point;
                min_neg_dist = std::abs(dist);
            }
        }
    }
    closest_pts.push_back(closest_neg);
    closest_pts.push_back(closest_pos);

    return closest_pts;
}



std::vector<Eigen::Vector2d> find_closest_points1(std::vector< Eigen::Vector2d > bdy_curve_pts, Eigen::Vector2d x0, Eigen::Vector2d unit_vec)
{
    double min_pos_dist = 100.;
    double min_neg_dist = 100.;

    std::vector<Eigen::Vector2d> closest_pts;
    std::vector<double> negative_distances;
    std::vector<Eigen::Vector2d> negative_points;
    std::vector<double> positive_distances;
    std::vector<Eigen::Vector2d> positive_points;

    #pragma omp parallel for
    for(int idx = 0; idx < bdy_curve_pts.size(); idx++)
    {

        Eigen::Vector2d point = bdy_curve_pts[idx];
        if( (point - x0).norm() > sqrt(2) )
        {
         continue;
        }
        double dist = signed_distance(point, x0, unit_vec);
        if(dist > 0.0)
        {
            positive_distances.push_back(dist);
            positive_points.push_back(point); 
        }
        else if (dist < 0.0)
        {
            negative_distances.push_back(dist);
            negative_points.push_back(point); 
        }

    }

    Eigen::Vector2d closest_positive_point, closest_negative_point;

    if (negative_distances.empty() && !positive_distances.empty())
    {
        auto closest_pos_it = std::min_element(positive_distances.begin(), positive_distances.end());
        int closest_pos_idx = std::distance(positive_distances.begin(), closest_pos_it);
        closest_positive_point = positive_points[closest_pos_idx];
        closest_pts.push_back(closest_positive_point);

        return closest_pts;
    }
    else if (positive_distances.empty() && !negative_distances.empty() )
    {
        auto closest_neg_it = std::max_element(negative_distances.begin(), negative_distances.end());
        int closest_neg_idx = std::distance(negative_distances.begin(), closest_neg_it);
        closest_negative_point = negative_points[closest_neg_idx];
        closest_pts.push_back(closest_negative_point);

        return closest_pts;
    }
    else if (!positive_distances.empty() && !negative_distances.empty() )
    {
        auto closest_pos_it = std::min_element(positive_distances.begin(), positive_distances.end());
        int closest_pos_idx = std::distance(positive_distances.begin(), closest_pos_it);
        closest_positive_point = positive_points[closest_pos_idx];
        closest_pts.push_back(closest_positive_point);

        auto closest_neg_it = std::max_element(negative_distances.begin(), negative_distances.end());
        int closest_neg_idx = std::distance(negative_distances.begin(), closest_neg_it);
        closest_negative_point = negative_points[closest_neg_idx];
        closest_pts.push_back(closest_negative_point);

        return closest_pts;
    }
    else
    {
        closest_positive_point = x0;
        closest_pts.push_back(closest_positive_point);
        return closest_pts;
    }


}


// End of the BoundaryNode class. Now we can move on to the other functions.



std::vector<std::vector<double>> make_grid(double x_min, double x_max, double y_min, double y_max, int num_pts_x, int num_pts_y) 
{

    std::vector<std::vector<double>> grid_pts;
    grid_pts.reserve(num_pts_x * num_pts_y);

    Eigen::VectorXd x_pts = Eigen::VectorXd::LinSpaced(num_pts_x, x_min, x_max);
    Eigen::VectorXd y_pts = Eigen::VectorXd::LinSpaced(num_pts_y, y_min, y_max);


    for (int i = 0; i < num_pts_x; ++i) 
    {
        for (int j = 0; j < num_pts_y ; ++j) 
        {
            grid_pts.push_back({x_pts[i], y_pts[j], 0.0});
        }
    }
    return grid_pts;

}


std::vector< Eigen::Vector2d > define_bdy_curve(int N_bdy_pts, double x_min, double x_max)
{
    double eps = 0.2;
    Eigen::VectorXd x_vals = Eigen::VectorXd::LinSpaced(N_bdy_pts, x_min - eps, x_max +eps);
    Eigen::VectorXd y_vals(N_bdy_pts);

    std::vector< Eigen::Vector2d > bdy_curve_pts;

    for (int i = 0; i < N_bdy_pts; i++)
    {
        y_vals[i] = ( 5*sin(M_PI/50 * x_vals[i]) + 10 );
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

    #pragma omp parallel for
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

        if (y_mesh_pt < y_curve) 
        {
            grid_pts[i][2] = 2;

            #pragma omp critical
            {
                solid_points.push_back({x_mesh_pt, y_mesh_pt});
            }
        }
    }

    return solid_points;
}


std::vector<BdyNode> label_bdy_points(std::vector<std::vector<double>>& grid_pts,
    const std::vector<std::vector<double>>& solid_points)
{
    std::vector<double> x_pts;
    std::vector<double> y_pts;
    std::vector<BdyNode> boundary_nodes;
    double eps = 1e-5;

    double dx = grid_pts[1][0] - grid_pts[0][0];
    double dy = grid_pts[1][1] - grid_pts[0][1];

    for (int i = 0; i < grid_pts.size(); i++)
    {
        std::vector< std::valarray<double> > velocity_vecs;
        if ( grid_pts[i][2] == 2)
        {
            continue; // If you encounter a solid point, move to next iteration of loop
        }
        
        //Eigen::VectorXd P(2); 
        //P << grid_pts[i][0], grid_pts[i][1];

        Eigen::Vector2d P(grid_pts[i][0], grid_pts[i][1]);

        for (int j = 0; j < solid_points.size(); j++)
        {
            //Eigen::VectorXd Q(2);
            //Q << solid_points[j][0], solid_points[j][1];

            Eigen::Vector2d Q(solid_points[j][0], solid_points[j][1]);

            Eigen::Vector2d vel_vec = Q - P;
            Eigen::Vector2d vel_vec2;


            if (vel_vec.norm() <=  sqrt(2.0)*std::max(dx, dy) + eps)
            {

                // velocity_vecs.push_back( (Q - P)/vector_norm(Q - P) );
                grid_pts[i][2] = 1;

                BdyNode x_b;
                x_b.Position = Eigen::Vector2d(P[0], P[1]);
                x_b.velocity_vecs.push_back( vel_vec.normalized() );
                std::vector<double> ez_read_vec(vel_vec.data(), vel_vec.data() + vel_vec.size());
                x_b.vecs.push_back(ez_read_vec);

                for(int k = 0; k < solid_points.size(); k++)
                {
                    if(k == j)
                    {continue;}
                    Eigen::Vector2d Q(solid_points[k][0], solid_points[k][1]);
                
                    vel_vec2 = Q - P;
                    if(vel_vec2.norm() <= std::sqrt(2)*std::max(dx, dy) + eps)
                    {
                        x_b.velocity_vecs.push_back( vel_vec2.normalized() );
                        Eigen::Vector2d V = vel_vec2.normalized();
                        std::vector<double> ez_read_vec(V.data(), V.data() + V.size());
                        x_b.vecs.push_back(ez_read_vec);
                    }

                }

                boundary_nodes.push_back(x_b);
                break;


            }


        } 

    }

    return boundary_nodes;

}

void distances_and_normals(std::vector<BdyNode>& boundary_nodes, std::vector< Eigen::Vector2d >& bdy_curve_pts)
{
    double eps = 3*1e-2;
    #pragma omp parallel for
    for(int j = 0; j < boundary_nodes.size(); j++)
    {
        Eigen::Vector2d x_b = boundary_nodes[j].Position;

        std::vector< Eigen::Vector2d > velocity_vecs = boundary_nodes[j].velocity_vecs;

        for(int i = 0; i < velocity_vecs.size(); i++)
        {
            Eigen::Vector2d unit_vel_vec = boundary_nodes[j].velocity_vecs[i];
            std::vector<Eigen::Vector2d> closest_pts = find_closest_points1(bdy_curve_pts, x_b, unit_vel_vec);
            Eigen::Vector2d pt1 = closest_pts[0];
            Eigen::Vector2d v = closest_pts[1] - closest_pts[0];
            v = v.normalized();

            Eigen::Matrix2d intersection_mat;
            intersection_mat << v[0], -unit_vel_vec[0],
                                v[1], -unit_vel_vec[1];

            Eigen::Vector2d intersection_rhs_vec;
            intersection_rhs_vec << x_b[0] - pt1[0], x_b[1] - pt1[1];

            Eigen::Vector2d intersection_distances = intersection_mat.colPivHouseholderQr().solve(intersection_rhs_vec);

            double delta = intersection_distances[1];

            if( delta > sqrt(2) + eps )
            {
                std::cout << "delta = " << delta << "\n\n";
                continue;
            }
            else
            {
                boundary_nodes[j].distances.push_back(delta);
                Eigen::Vector2d normal_vec(1, -v[0]/v[1]);
                normal_vec = normal_vec/normal_vec.norm();
                boundary_nodes[j].normals.push_back(normal_vec);

                
            }

        }

    }
}

void constructMatrix(std::vector<BdyNode>& boundary_nodes)
{

}

void write_data(const std::vector<std::vector<double>>& grid_pts, const std::vector< Eigen::Vector2d >& bdy_curve_pts,
    const std::vector<std::vector<double>>& solid_points, std::vector<BdyNode>& boundary_nodes )
{
    const int width = 10;
    const int precision = 5;

    FILE* outfile1 = fopen("boundary_pts.dat", "w");
    std::vector<std::vector<double>> boundary_node_positions;
    for(int i = 0; i < boundary_nodes.size(); i++)
    {
        Eigen::Vector2d x_b = boundary_nodes[i].Position;
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


    double alpha = 0.2;
    double R = 1;
    int N_bdy_pts = 3200;
    int N_repeats = 1;
    int n_grid_pts = 100;
    double x_min= 0;
    double x_max = 200;
    double y_min = 0;
    double y_max = 50;
    int num_pts_x = 200;
    int num_pts_y = 50;


    std::vector<std::vector<double>> grid_pts = make_grid(x_min, x_max, y_min, y_max, num_pts_x, num_pts_y); // Won't need this, grid points comes from LBM code.

    std::vector< Eigen::Vector2d > bdy_curve_pts = define_bdy_curve(N_bdy_pts, x_min, x_max);

    std::vector<std::vector<double>> solid_points = label_solid_pts(grid_pts, bdy_curve_pts); // Will need this. 

    std::vector<BdyNode> boundary_nodes = label_bdy_points(grid_pts, solid_points);

    distances_and_normals(boundary_nodes, bdy_curve_pts);

    for( BdyNode node : boundary_nodes )
    {
        int n = node.velocity_vecs.size();
        int m = node.normals.size();
        if(n != m)
        printf("Houston, we have a problem \n\n");
    }

    write_data(grid_pts, bdy_curve_pts, solid_points, boundary_nodes );

    constructMatrix(boundary_nodes);


}