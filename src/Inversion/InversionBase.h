#ifndef INVERSION_BASE_H
#define INVERSION_BASE_H
#include <cassert>
#include <functional>
#include <random>

#include "Fwd.h"
#include "gs.h"
#include "linterp.h"
using namespace std;
using namespace Eigen;
class InversionBase : public Fwd {
   public:
    InversionBase();
    InversionBase(const Mesh& mesh_, const Observation& ob_,
                  unsigned long long field_flag_ = Compute_g_r);
    ~InversionBase();

    virtual void invert() = 0;

    void output_obs_data(string out_name);

    void output_predicted_data(string out_name);

    VectorXd& get_result() { return this->m; }
    VectorXd get_predicted_field();
    double get_final_misfit() const { return final_misfit; }
    double get_final_lambda() const { return final_lambda; }

    void set_Wd(const VectorXd& sigma);

    void set_dobs(const VectorXd& d);
    // void set_dobs(const VectorXd &d, double relative_error);
    // //relative_error, e.g. 0.05
    void set_dobs2(const VectorXd& d,
                   double constant_error);  // relative_error, e.g. 0.05
    void set_dobs(const VectorXd& d, double relative_error,
                  double constant_error = 0);
    void set_dobs(const VectorXd& d, vector<double> relative_error,
                  vector<double> constant_error);

    void set_depth_weighting(double beta = 2, int flag = 0);
    void set_target_misfit(double x) { this->target_misfit = x; }
    void set_lambda_decreasing_rate(double rate) {
        this->lambda_decreasing_rate = rate;
    }
    void set_max_lambda(double max_lambda_) { this->max_lambda = max_lambda_; }
    void set_n_lambda(int n_lambda_) { this->n_lambda = n_lambda_; }

    void show_differece_matrix(unsigned int direction);

    // void set_weights_of_objectives(double as, double ar, double a_theta,
    // double a_phi);
    void set_weights_of_objectives(double as, double ar, double a_theta,
                                   double a_phi, double a_crg = 0);

    // set mesh of tesseroids
    void set_mesh(const Mesh& mesh0);

    // set computation points, overwrite
    void set_observation(const Observation& ob0);

    void set_density_to_mesh();
    void set_reference_model_to_mesh();
    void set_min_max_to_mesh();

    void set_reference_model(VectorXd& m_ref);
    void set_geometry_reference_model(VectorXd& m_ref);
    void set_m(VectorXd& m_);
    void set_m_ini(VectorXd& m_ini_);
    void set_min_max(VectorXd& m_min_, VectorXd& m_max_);

    void set_petrophysics_constraint(VectorXd& s,
                                     function<double(double)> relation);

    /**
     * @brief  set the 3D interpolater of the structural constraint model using
     * gridded data from a text file. File format requirement:
     * The first 3 colums specify coordinates and the fourth column specifies
     * model values. Latitude ranges from -90 (south) to 90 (north), lontigude
     * ranges from -180 to 180, depth is greater than 0, with positive axis
     * pointing downwards.
     *
     * @param  filename  a file containing a 3D a-priori model
     * @param  lat_size  number of grid nodes along latitude dimension
     * @param  lon_size  number of grid nodes along longitude dimension
     * @param  dep_size  number of grid nodes along depth dimension
     * @param  format_of_coordinates  any combination of characters 'x', 'y' and
     * 'z', specifying how the coordinates are ordered. "yxz" means:
     * longitude(column 1) latitude(column 2)  depth(column 3)  value(column 4)
     * "xyz" means:
     * latitude(column 1) longitude(column 2) depth(column 3) value(column 4)
     * "zxy" means:
     * depth(column 1)  latitude(column 2) longitude(column 3) value(column 4)
     * "zyx" means:
     * depth(column 1)  longitude(column 2) latitude(column 3) value(column 4)
     * "xzy" means:
     * latitude(column 1)  depth(column 2) longitude(column 3) value(column 4)
     * "yzx" means:
     * longitude(column 1) depth(column 2) latitude(column 3) value(column 4)
     *
     * @param  fast_dimension  0 or 1, specifying which dimension changes
     * fastest. 0 represents latitude, 1 represents longitude. Whatever
     * fast_dimension is, the depth dimension is assumed to be changing slowest.
     */
    void create_crg_model_from_data(string filename, int lat_size, int lon_size,
                                    int dep_size,
                                    string format_of_coordinates = "yxz",
                                    int fast_dimension = 0,
                                    double reference_level = 6378137);

    void create_ref_model_from_data(string filename, int lat_size, int lon_size,
                                    int dep_size,
                                    string format_of_coordinates = "yxz",
                                    int fast_dimension = 0,
                                    double reference_level = 6378137);

    void set_interpolator_m0s(InterpMultilinear<3, double>* interp) {
        this->use_cross_gradient_constraint = true;
        interpolator_m0s = interp;
    }

    void set_interpolator_m0(InterpMultilinear<3, double>* interp) {
        this->use_petrophysical_constraint = true;
        interpolator_m0 = interp;
    }

    void set_constraint_region(double depth[2], double lat[2], double lon[2]) {
        for (int i = 0; i < 2; i++) {
            this->constraint_depth[i] = depth[i];
            this->constraint_lat[i] = lat[i];
            this->constraint_lon[i] = lon[i];
        }
    }

    void cumulative_sensitivity();
    void result2vtk(string filename);
    void result2txt(string filename);
    void result2netcdf(string filename);

    int number_unknowns() { return Nm; }

   protected:
    void init_matrices();

    void set_S();
    void set_Tmatrix();
    void set_Vmatrix();

    void set_difference_matrix();

    void update_S_crg();

    void out_data(const VectorXd& d, string out_name);

   protected:
    VectorXd m;
    VectorXd dobs;
    VectorXd m0;    // reference model
    VectorXd m0_s;  // structure similarity
    VectorXd m_min;
    VectorXd m_max;

    // a 3D interpolator of cross-gradient constraint model. Once the mesh is
    // changed, the constraint model is updated using this interpolator
    InterpMultilinear<3, double>* interpolator_m0;

    // a 3D interpolator of the reference model
    InterpMultilinear<3, double>* interpolator_m0s;

    bool use_cross_gradient_constraint;
    bool use_petrophysical_constraint;

    // extent of the constraint region
    double constraint_depth[2];
    double constraint_lat[2];
    double constraint_lon[2];

    SMatrix Wd;   // diagonal matrix whose i-th element is 1/sigma_i, sigma_i is
                  // the standard deviation of i-th datum
    SMatrix R_r;  // radius weighting
    SMatrix V;

    SMatrix S_s;
    SMatrix S_r;
    SMatrix S_theta;
    SMatrix
        S_phi;  // spatially dependent 3D weighting function, diagonal matrix
    SMatrix S_crg;

    SMatrix D_s;
    SMatrix D_r1;
    SMatrix D_theta1;
    SMatrix D_phi1;  // first-order finite-difference representation

    SMatrix T_r;
    SMatrix T_theta;
    SMatrix T_phi;
    double a_s, a_r, a_theta, a_phi, a_crg;

    double final_misfit;
    double final_lambda;

    double depth_weighting_factor;
    double target_misfit;
    double max_lambda;
    double lambda_decreasing_rate;
    int n_lambda;
    double reference_surface;

    VectorXd m_ini;
};
#endif
