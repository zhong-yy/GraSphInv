#include <chrono>
#include <fstream>
#include <random>
// #include<function>

#include "GaussNewtonInversion.h"
#include "timer.h"
void out_field2(ofstream& os2, Observation& ob, VectorXd& d_obs) {
  for (int i = 0; i < ob.get_n_obs(); i++) {
    const Point& p = ob(i);
    os2 << fixed;
    os2 << setw(30) << setprecision(7) << left << p(0) << setw(30) << left
        << 90.0 - p(1) * 180.0 / GS::PI << setw(30) << left
        << p(2) * 180.0 / GS::PI - 180.0;
    os2 << scientific;
    os2 << setw(30) << setprecision(15) << left << d_obs(i) << endl;
  }
}
int main() {
  double lats[2] = {10, 40};
  double lons[2] = {90, 120};
  double dep[2] = {0, 500000};
  double reference_surface = 6378137;

  Mesh tmesh;
  tmesh.generate_regular_geographic_mesh(lats, 15, lons, 15, dep, 5,
                                         reference_surface);

  double block_dep[2] = {100000, 300000};
  double block_lat[2] = {20, 30};
  double block_lon[2] = {100, 110};
  tmesh.set_parameter_in_a_region(block_lat[0], block_lat[1], block_lon[0],
                                  block_lon[1], block_dep[0], block_dep[1], 500,
                                  reference_surface);

  // tmesh.set_block_density(12, 27, 12, 27, 6, 15, 500);//
  // tmesh.set_block_density(6, 13, 6, 13, 6, 15, 500); //

  tmesh.out_model_vtk("test_model.vtk");
  tmesh.out_model_netcdf("test_model.nc");

  /******************Foward modelling*****************/
  VectorXd rho;
  tmesh.get_model_parameter_from_mesh(rho, 0);
  Observation ob;

  int nlat = 31;
  int nlon = 31;
  double observation_r = 6408137;
  // ob.generate(1, observation_r, observation_r, nlat,
  //             (theta_model[0]) * GS::PI / 180.0,
  //             (theta_model[1]) * GS::PI / 180.0, nlon,
  //             (phi_model[0]) * GS::PI / 180.0, (phi_model[1]) * GS::PI /
  //             180.0);
  ob.generate_geographic(1, observation_r, observation_r, nlat, lats[0],
                         lats[1], nlon, lons[0], lons[1]);
  ofstream os("sites");
  os << ob;

  Fwd forward(tmesh, ob,
              Compute_g_r |Compute_g_phi|Compute_g_theta| Compute_T_rr | Compute_T_rtheta | Compute_T_rphi |
                  Compute_T_thetatheta | Compute_T_thetaphi | Compute_T_phiphi,
              16);
  // Fwd forward(&tmesh, &ob, Compute_T_rr|Compute_T_rtheta, 8);

  Timer timer;
  timer.start();
  cout << "Generating synthetic data using volume integrals..." << endl;
  forward.set_integral_kernel_type(1);
  forward.compute_G();
  const Eigen::MatrixXd& G = forward.get_G();
  Eigen::VectorXd d_obs = G * rho;
  cout << "Calculation of synthetic data completed" << endl;
  timer.stop();
  cout << "Time: " << timer.getElapsedTimeInSec() << " s" << endl;
  GaussNewtonInversion inv(tmesh, ob,
                           Compute_g_r|Compute_g_phi|Compute_g_theta | Compute_T_rr | Compute_T_rtheta |
                               Compute_T_rphi | Compute_T_thetatheta |
                               Compute_T_thetaphi | Compute_T_phiphi);
  inv.set_dobs(d_obs);
  inv.output_obs_data("dobs_v");

  // surface integral
  Fwd forward2(tmesh, ob,
              Compute_g_r |Compute_g_phi|Compute_g_theta| Compute_T_rr | Compute_T_rtheta | Compute_T_rphi |
                  Compute_T_thetatheta | Compute_T_thetaphi | Compute_T_phiphi,
              16);
  timer.start();
  cout << "Generating synthetic data using surface integrals..." << endl;
  forward2.set_integral_kernel_type(0);
  forward2.compute_G();
  const Eigen::MatrixXd& Gs = forward2.get_G();
  Eigen::VectorXd d_obs_surf = Gs * rho;
  cout << "Calculation of synthetic data completed" << endl;
  timer.stop();
  cout << "Time: " << timer.getElapsedTimeInSec() << " s" << endl;
  GaussNewtonInversion inv2(tmesh, ob,
                            Compute_g_r |Compute_g_phi|Compute_g_theta| Compute_T_rr | Compute_T_rtheta |
                                Compute_T_rphi | Compute_T_thetatheta |
                                Compute_T_thetaphi | Compute_T_phiphi);
  inv2.set_dobs(d_obs_surf);
  inv2.output_obs_data("dobs_s");

  return 0;
}
