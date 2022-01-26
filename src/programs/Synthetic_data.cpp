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
    double lats[2] = {20, 40};
    double lons[2] = {90, 110};
    double dep[2] = {0, 400000};
    double reference_surface = 6378137;

    Mesh tmesh;
    tmesh.generate_regular_geographic_mesh(lats, 40, lons, 40, dep, 20,
                                           reference_surface);

    double block_dep[2] = {80000, 200000};
    double block_lat[2] = {27, 33};
    double block_lon[2] = {95, 98};
    tmesh.set_parameter_in_a_region(block_lat[0], block_lat[1], block_lon[0],
                                    block_lon[1], block_dep[0], block_dep[1],
                                    200, reference_surface);

    double block_dep2[2] = {80000, 200000};
    double block_lat2[2] = {32, 35};
    double block_lon2[2] = {102, 105};
    tmesh.set_parameter_in_a_region(block_lat2[0], block_lat2[1], block_lon2[0],
                                    block_lon2[1], block_dep2[0], block_dep2[1],
                                    200, reference_surface);

    double block_dep3[2] = {80000, 200000};
    double block_lat3[2] = {25, 28};
    double block_lon3[2] = {102, 105};
    tmesh.set_parameter_in_a_region(block_lat3[0], block_lat3[1], block_lon3[0],
                                    block_lon3[1], block_dep3[0], block_dep3[1],
                                    -200, reference_surface);

    // tmesh.set_block_density(12, 27, 12, 27, 6, 15, 500);//
    // tmesh.set_block_density(6, 13, 6, 13, 6, 15, 500); //

    tmesh.out_model_vtk("test_model.vtk");
    tmesh.out_model_netcdf("test_model.nc");

    Mesh tmesh2;
    tmesh2.generate_regular_geographic_mesh(lats, 40, lons, 40, dep, 20,
                                            reference_surface);
    tmesh2.set_parameter_in_a_region(block_lat[0], block_lat[1], block_lon[0],
                                     block_lon[1], block_dep[0], block_dep[1],
                                     15, reference_surface);
    ofstream out_cross_gradient_constraint("crg_model");
    for (int k = tmesh2.get_nr_level_0() - 1; k >= 0; k--) {
        for (int j = 0; j < tmesh2.get_nphi_level_0(); j++) {
            for (int i = tmesh2.get_ntheta_level_0() - 1; i >= 0; i--) {
                Cell* c = tmesh2.get_element_level_0(i, j, k);
                double rc, thetac, phic;
                c->get_center(rc, thetac, phic);
                double latc = 90.0 - thetac * 180.0 / GS::PI;
                double lonc = phic * 180 / GS::PI - 180;
                out_cross_gradient_constraint
                    << setw(15) << fixed << setprecision(7) << left << latc
                    << setw(15) << setprecision(7) << left << lonc << setw(15)
                    << setprecision(3) << left
                    << (reference_surface - rc) / 1000 << scientific
                    << setprecision(15) << left << c->get_parameter() << endl;
            }
        }
    }

    // function<double(double)> relation = [](double x) -> double {
    //   return (2.0 * x + 50.0);
    // };  //
    // ofstream out_petrophysical_constraint("ref_model");
    // for (int k = tmesh2.get_nr_level_0() - 1; k >= 0; k--) {
    //   for (int j = 0; j < tmesh2.get_nphi_level_0(); j++) {
    //     for (int i = tmesh2.get_ntheta_level_0() - 1; i >= 0; i--) {
    //       Cell* c = tmesh2.get_element_level_0(i, j, k);
    //       double rc, thetac, phic;
    //       c->get_center(rc, thetac, phic);
    //       double latc = 90.0 - thetac * 180.0 / GS::PI;
    //       double lonc = phic * 180.0 / GS::PI - 180.0;
    //       out_petrophysical_constraint
    //           << setw(15) << fixed << setprecision(7) << left << latc <<
    //           setw(15)
    //           << setprecision(7) << left << lonc << setw(15) <<
    //           setprecision(3)
    //           << left << (reference_surface - rc) / 1000.0 << scientific
    //           << setprecision(15) << left << relation(c->get_parameter()) <<
    //           endl;
    //     }
    //   }
    // }

    /******************Foward modelling*****************/
    VectorXd rho;
    tmesh.get_model_parameter_from_mesh(rho, 0);
    Observation ob;

    int nlat = 41;
    int nlon = 41;
    double observation_r = 6428137;
    // ob.generate(1, observation_r, observation_r, nlat,
    //             (theta_model[0]) * GS::PI / 180.0,
    //             (theta_model[1]) * GS::PI / 180.0, nlon,
    //             (phi_model[0]) * GS::PI / 180.0, (phi_model[1]) * GS::PI /
    //             180.0);
    ob.generate_geographic(1, observation_r, observation_r, nlat, lats[0],
                           lats[1], nlon, lons[0], lons[1]);
    ofstream os("sites");
    os << ob;

    Fwd forward(tmesh, ob, Compute_g_r, 4);
    // Fwd forward(tmesh, ob,
    // Compute_T_rr|Compute_T_rtheta|Compute_T_rphi|Compute_T_thetatheta|Compute_T_thetaphi|Compute_T_phiphi,8);
    // Fwd forward(&tmesh, &ob, Compute_T_rr|Compute_T_rtheta, 8);

    Timer timer;
    timer.start();
    cout << "Generating synthetic data..." << endl;
    forward.compute_G();
    const Eigen::MatrixXd& G = forward.get_G();
    Eigen::VectorXd d_obs = G * rho;
    cout << "Calculation of synthetic data completed" << endl;
    timer.stop();
    cout << "Time: " << timer.getElapsedTimeInSec() << " s" << endl;

    VectorXd noise;
    noise.resize(d_obs.rows());
    double equipment_noise = 0.0;
    static normal_distribution<double> normal_dist(0, 1);
    static default_random_engine e(time(0));
    //  double std=0.01*d_obs.maxCoeff();
    VectorXd abs_dobs = d_obs.cwiseAbs();

    cout << "||d_obs||=" << d_obs.norm() << endl;
    cout << "0.01*max(d_obs)=" << 0.01 * abs_dobs.maxCoeff() << endl;

    equipment_noise = 0.01 * abs_dobs.maxCoeff();

    for (int i = 0; i < d_obs.rows(); i++) {
        //    noise(i) = (0.02 * fabs(d_obs(i)) + equipment_noise) *
        //    normal_dist(e);
        noise(i) = (0.02 * fabs(d_obs(i)) + equipment_noise) * normal_dist(e);
        //    noise(i)=std*normal_dist(e);
        d_obs(i) = d_obs(i) + noise(i);
    }

    // int n_fields=forward.get_n_fields();
    // int n_obs=forward.get_N_obs();
    // cout<<n_fields<<endl;
    // cout<<n_obs<<endl;
    // for(int i=0;i<n_fields;i++){
    //   VectorXd temp=d_obs.segment(i*n_obs,n_obs);
    //   VectorXd abs_dobs2=temp.cwiseAbs();
    //   equipment_noise=0.01*abs_dobs2.maxCoeff();
    //   cout<<"0.01*max(d_obs)="<<0.01*abs_dobs2.maxCoeff()<<endl;
    //   for(int j=0;j<d_obs.rows();j++){
    //     noise(i) = (0.02 * fabs(d_obs(i)) + equipment_noise) *
    //     normal_dist(e);
    // //    noise(i)=std*normal_dist(e);
    //     d_obs(i) = d_obs(i) + noise(i);
    //   }
    // }
    // for()

    // cout << d_obs.rows() << endl;

    // GaussNewtonInversion inv(tmesh, ob,
    // Compute_T_rr|Compute_T_rtheta|Compute_T_rphi|Compute_T_thetatheta|Compute_T_thetaphi|Compute_T_phiphi);
    GaussNewtonInversion inv(tmesh, ob, Compute_g_r);

    inv.set_dobs(d_obs);
    inv.output_obs_data("dobs");

    return 0;
}
