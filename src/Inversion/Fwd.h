#ifndef _FWD
#define _FWD
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>

#include <algorithm>
#include <bitset>

#include "GravityField.h"
#include "Mesh.h"
#include "Observation.h"
#include "gs.h"
using namespace std;

class Fwd {
   public:
    Fwd();
    Fwd(const Mesh &mesh_, const Observation &ob_,
        unsigned long long field_flag_ = 0x03ff, int GLQ_order_ = 2);
    Fwd(const Mesh &mesh_, const Observation &ob_, bitset<32> field_flag,
        int GLQ_order_ = 2);
    virtual ~Fwd();

    VectorXd compute_gobs(const VectorXd &rho) {
        VectorXd d;
        if (use_wavelet) {
            this->G_vec_mul(rho, d);
        } else {
            d = this->G * rho;
        }
        return d;
    }

    void set_field_flag(unsigned long long field_flag1) {
        field_flag = field_flag1;
    }

    void GT_vec_mul(const VectorXd &vec, VectorXd &product) const;
    void G_vec_mul(const VectorXd &vec, VectorXd &product) const;

    void compute_G();
    void compute_G_wavelet();

    // bind mesh of tesseroids
    void set_mesh(const Mesh &mesh0);

    // bind computation points
    void set_observation(const Observation &ob0);

    void set_use_wavelet(bool use_wavelet0);

    const Eigen::MatrixXd &get_G() { return G; }

    unsigned int get_n_fields() const { return field_flag.count(); }

    unsigned int get_N_obs() const { return N_obs; }

    void set_GLQ_order(int order) { this->GLQ_order = order; }

    void set_integral_kernel_type(int ikt) { integral_kernel_type = ikt; }

    Mesh &get_mesh() { return mesh; }
    void set_compression_threshold(double eps) {
        assert((std::fabs(eps) < 1e-15 || eps > 0) && (eps < 1));
        this->compression_threshold = eps;
    }

   protected:
    Eigen::MatrixXd G;
    Mesh mesh;
    Observation ob;

    bool use_wavelet;
    int Nm;     // number of parameters
    int N_obs;  // number of data points
    int Nd;     // number of observations

    int n_dy_for_wavelet;

    int GLQ_order;
    int integral_kernel_type;

    bitset<32> field_flag;

    // relative threshold for wavelet compression, default 0.005
    double compression_threshold;

    vector<vector<int>> comp_col_ids;
    vector<vector<double>> comp_G_coeffs;

    int inverse_wavelet_transform_vec(vector<double> &data) const;
    int compress_vec(vector<double> &data, vector<int> &ids,
                     vector<double> &coeffs, double relative_threshold);
    int wavelet_transform_vec(vector<double> &data) const;

    void display_info_fields() const;
    /**
     * 0 V
     * 1 gr
     * 2 g_theta
     * 3 g_phi
     * 4 T_rr
     * 5 T_rtheta
     * 6 T_rphi
     * 7 T_thetatheta
     * 8 T_thetaphi
     * 9 T_phiphi
     */
};
#endif
