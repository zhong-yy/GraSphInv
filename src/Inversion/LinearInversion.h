#ifndef LINEAR_INVERSION_H
#define LINEAR_INVERSION_H

#include "InversionBase.h"
class LinearInversion : public InversionBase {
public:
  LinearInversion(){};
  LinearInversion(const Mesh& mesh_,
                  const Observation& ob_,
                  unsigned long long field_flag_)
      : InversionBase(mesh_, ob_, field_flag_) {}
  void invert() override;

  void set_CG_parameter(double cg_tol_, double cg_iteration_factor_) {
    cg_tol = cg_tol_;
    cg_iteration_factor_ = abs(cg_iteration_factor_);
    this->cg_iteration_factor = cg_iteration_factor_;
  }

 protected:
  double cg_tol;               // tolerance for conjugate gradient
  double cg_iteration_factor;  // CG iteration number
  double stag_tol;             // stagnation
};
#endif
