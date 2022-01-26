#ifndef _MESH
#define _MESH

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <netcdf>  //The library netcdf-cxx is required
#include <set>
#include <unordered_map>
#include <vector>

#include "Cell.h"
using namespace netCDF;
using namespace netCDF::exceptions;
#define NC_ERR 2

class Mesh {
   public:
    Mesh();
    Mesh(const Mesh&);
    Mesh& operator=(const Mesh&);
    ~Mesh();

    void clear_all();

    void show_ordering();
    void set_n_parameter(int n) {
        this->n_parameters = n;
        for (int i = 0; i < leaf_cells.size(); i++) {
            leaf_cells[i]->parameters.resize(n);
            for (int j = 0; j < n; j++) {
                leaf_cells[i]->set_parameter(0, j);
            }
        }
    }

    /**
     * @brief To generate a regular mesh, similar to generate_regular_mesh(...)
     *        but geographic coordinates are used here
     *
     * @param lats       Range of latitude
     * @param n_lat      How many parts the mesh will be divided from south to
     * north
     * @param lons       Range of longitude
     * @param n_lon      How many parts the mesh will be divided from west to
     * east
     * @param depth_model Range of depth
     * @param n_depth    How many parts the mesh will be divided in depth
     * direction
     * @param ref_suface Reference surface
     * @param num_para How many kinds of model parameters will be stored
     */
    void generate_regular_geographic_mesh(double lats[2], int n_lat,
                                          double lons[2], int n_lon,
                                          double depth_model[2], int n_depth,
                                          double ref_suface = 6378137,
                                          int num_para = 1);

    /**
     * @brief To generate a regular mesh
     *
     * @param theta_model  Range of theta, in degree. Domain: [0,180]
     * @param n_theta      How many parts the mesh will be divided in theta axis
     * @param phi_model    Range of phi, in degree. Domain: [0, 360]
     * @param n_phi        How many parts the mesh will be divided in phi axis
     * @param r_model      Range of radius, in meter
     * @param n_r          How many parts the mesh will be divided in radial
     * axis
     * @param num_para     How many kinds of model parameters will be stored
     *
     * @note theta is 0deg at the north pole, 180deg at the south pole, and phi
     * is 0deg at prime meridian
     */
    void generate_regular_mesh(double theta_model[2], int n_theta,
                               double phi_model[2], int n_phi,
                               double r_model[2], int n_r, int num_para = 1);

    void generate_regular_mesh(double theta_model[2], int n_theta,
                               double phi_model[2], int n_phi,
                               VectorXd& r_points, int num_para = 1);

    //层厚度按等差数列方式递增
    void generate_radially_varying_mesh(double theta_model[2], int n_theta,
                                        double phi_model[2], int n_phi,
                                        double h1, double depth, int n_r,
                                        double reference_surface = 6378137,
                                        int num_para = 1);

    void generate_radially_varying_mesh_pow(double theta_model[2], int n_theta,
                                            double phi_model[2], int n_phi,
                                            double depth_1, double depth_2,
                                            int n_r, double nth_power = 2,
                                            double reference_surface = 6378137,
                                            int num_para = 1);

    void generate_radially_varying_mesh_log(double theta_model[2], int n_theta,
                                            double phi_model[2], int n_phi,
                                            double depth_1, double depth_2,
                                            int n_r, double base = 2,
                                            double reference_surface = 6378137,
                                            int num_para = 1);

    void set_parameter_in_a_region(double lat0, double lat1, double lon0,
                                   double lon1, double d0, double d1,
                                   double para_value,
                                   double reference_surface = 6378137,
                                   int i_th = 0);

    void set_block_parameter(unsigned int i_min, unsigned int i_max,
                             unsigned int j_min, unsigned int j_max,
                             unsigned int k_min, unsigned int k_max,
                             double para_value, int i_th = 0);

    void set_ith_cell_parameter(double para_value, int id, int i) {
        leaf_cells[id]->set_parameter(para_value, i);
    };

    // return a key-value map, whose key is the id of the cell that is to be
    // subdivided, and the value is the pointer to the cell
    map<unsigned int, Cell*> refinement(int i);
    map<unsigned int, Cell*> refinement(Cell* c);

    void refine_cell_across_a_interface(Cell* c, double r_interface,
                                        double r_resolution, double upper_part,
                                        double lower_part, int i_para = 0);

    void refine_cell_across_2_interfaces(Cell* c, double r_interface1,
                                         double r_interface2,
                                         double r_resolution, double upper_part,
                                         double middle_part, double lower_part,
                                         int i_para = 0);
    void out_model_txt(string filename, int ith_para = 0,
                       double reference_surface = 6378137);
    void out_model_vtk(
        string filename, int n = 1,
        vector<string> parameter_name = vector<string>(1, "density"));
    void out_model_vtk_linear_projection(
        string filename, int n = 1, double reference_surface = 6378137,
        vector<string> parameter_name = vector<string>(1, "density"));
    void convert_txt_2_vtk(string txtfile, string vtkfile, int n);
    int out_model_netcdf(string filename, int ith_para = 0,
                         string VAL_NAME = "density",
                         string VAL_UNITS = "kg/m3");
    void fill_data(int offset_i, int offset_j, int offset_k, double*** data,
                   Cell* c, int max_level, int ith_para);

    int n_elems() const { return leaf_cells.size(); }

    void rearrange_id();

    bool great_equal(long double left, long double right);

    void rthetaphi_xyz(long double r, long double theta, long double phi,
                       long double& x, long double& y, long double& z);

    void get_model_parameter_from_mesh(VectorXd& m, int ith = 0);

    void sort(int level);

    int get_reordered_id(const unsigned int& i) const;
    Tesseroid& get_elem(unsigned int i);
    Tesseroid& get_elem_reordered(unsigned int i);

    Cell* get_element_level_0(unsigned int i_lat, unsigned int j_lon,
                              unsigned int k_r);
    int get_nr_level_0() const { return nr; }
    int get_ntheta_level_0() const { return ntheta; }
    int get_nphi_level_0() const { return nphi; }
    //   protected:
    // vector<Cell *> root_cells;
    // vector<Face *> root_faces;

    void get_minimum_size(double& dlat, double& dlon, double& dr, int& lev) {
        assert(cells[cells.size() - 1].size() > 0);
        Cell* c = cells[cells.size() - 1][0];
        lev = c->get_level();
        c->get_size(dr, dlat, dlon);
        dlon *= 180. / GS::PI;
        dlat *= 180. / GS::PI;
    }
    int n_parameters;

    // number of cells along x, y, z directions for level 0
    int nr, ntheta, nphi;

    int num_leaf_cells;
    int num_leaf_faces;

    double r_lim[2], theta_lim[2], phi_lim[2];

    vector<Cell*> leaf_cells;
    vector<Face*> leaf_faces;
    vector<Cell*> ro_leaf_cells;

    vector<vector<Cell*>> cells;  // cells of different levels
    vector<vector<Face*>> faces;  // faces of different levels

    vector<int> num_cell;
    vector<int> num_face;

    // ordering model elements in a continuous way (to make consecutive indices
    // label adjacent model cells)
    vector<int> new_to_old_index;
    vector<int> old_to_new_index;
    // cells which which are re-ordered in a continuous way

    // vector<double> r_intervals;
    VectorXd r_points;
};
#endif
