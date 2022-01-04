#ifndef GN_INVERSION_WAVELET_H
#define GN_INVERSION_WAVELET_H
#include "AdaptiveInversion.h"
class GN_Inversion_Wavelet : public AdaptiveInversion {
   public:
    GN_Inversion_Wavelet()
        : AdaptiveInversion(),
          record_process(false),
          cg_tol(1e-6),
          stag_tol(0.025),
          GN_iter(10) {}
    GN_Inversion_Wavelet(const Mesh& mesh_, const Observation& ob_,
                         unsigned long long field_flag_)
        : AdaptiveInversion(mesh_, ob_, field_flag_),
          record_process(false),
          cg_tol(1e-6),
          stag_tol(0.025),
          GN_iter(10) {}

    void invert() override;

    void record_every_iteration() { this->record_process = true; }

    // when cg_iteration_factor is greater than 1, it's the iteration
    // number of conjugate gradient method; when it's lower than 1, CG
    // iteration number is  cg_iteration_factor times the cells number
    void set_CG_parameter(double tol, double cg_iteration_factor);
    void set_stagnation_tolerance(double x) { stag_tol = x; }
    void set_max_GN_iterations(int x) { GN_iter = x; }

    void display_inversion_parameters() const;

    void get_curvature(const vector<double>& x, const vector<double>& y,
                       vector<double>& curvature);

    VectorXd solve_cg(const double tol, const int& maxit, double& error,
                      int& iterations, const double& lambda, const VectorXd& mk,
                      const VectorXd& mref, const SMatrix& Wd,
                      const VectorXd& dobs, const SMatrix& P, const SMatrix& Ws,
                      const SMatrix& Wtheta, const SMatrix& Wphi,
                      const SMatrix& Wr, const SMatrix& Ttheta,
                      const SMatrix& Tphi, const SMatrix& Tr);

   protected:
    bool record_process;
    double cg_tol;               // tolerance for conjugate gradient
    double cg_iteration_factor;  // CG iteration number
    double stag_tol;             // stagnation
    int GN_iter;                 // Gauss-Newton iteration number
};

#endif
