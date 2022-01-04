#ifndef GAUSSNEWTONINVERSION_H
#define GAUSSNEWTONINVERSION_H
#include "AdaptiveInversion.h"
class GaussNewtonInversion : public AdaptiveInversion {
   public:
    GaussNewtonInversion()
        : AdaptiveInversion(),
          record_process(false),
          cg_tol(1e-6),
          stag_tol(0.025),
          GN_iter(10),
          method_id(1) {}
    GaussNewtonInversion(const Mesh& mesh_, const Observation& ob_,
                         unsigned long long field_flag_)
        : AdaptiveInversion(mesh_, ob_, field_flag_),
          record_process(false),
          cg_tol(1e-6),
          stag_tol(0.025),
          GN_iter(10),
          method_id(1) {}

    void invert() override;

    void invert_with_Eigen_CG();
    void invert_with_own_CG();

    void record_every_iteration() { this->record_process = true; }

    void set_method_id(int method_id0) { this->method_id = method_id0; }
    // when cg_iteration_factor is greater than 1, it's the iteration
    // number of conjugate gradient method; when it's lower than 1, CG
    // iteration number is  cg_iteration_factor times the cells number
    void set_CG_parameter(double tol, double cg_iteration_factor);
    void set_stagnation_tolerance(double x) { stag_tol = x; }
    void set_max_GN_iterations(int x) { GN_iter = x; }

    VectorXd solve_cg(const double tol, const int& maxit, double& error,
                      int& iterations, const double& lambda, const VectorXd& mk,
                      const VectorXd& mref, const SMatrix& Wd,
                      const VectorXd& dobs, const SMatrix& P, const SMatrix& Ws,
                      const SMatrix& Wtheta, const SMatrix& Wphi,
                      const SMatrix& Wr, const SMatrix& Ttheta,
                      const SMatrix& Tphi, const SMatrix& Tr);

    // The G matrix is not stored explicitly. Every element of G is calculated
    // every time G*m or G^T*x is required.
    // VectorXd solve_cg_without_G_matrix(const double tol, const int& maxit,
    //                                   double& error, int& iterations,
    //                                   const double& lambda, const VectorXd&
    //                                   mk, const VectorXd& mref, const
    //                                   SMatrix& Wd, const VectorXd& dobs,
    //                                   const SMatrix& P, const SMatrix& Ws,
    //                                   const SMatrix& Wtheta, const SMatrix&
    //                                   Wphi, const SMatrix& Wr, const SMatrix&
    //                                   Ttheta, const SMatrix& Tphi, const
    //                                   SMatrix& Tr);

    void display_inversion_parameters() const;

    void get_curvature(const vector<double>& x, const vector<double>& y,
                       vector<double>& curvature);

   protected:
    bool record_process;
    int GN_iter;  // Gauss-Newton iteration number
    int method_id;
    double cg_tol;               // tolerance for conjugate gradient
    double cg_iteration_factor;  // CG iteration number
    double stag_tol;             // stagnation
};

#endif
