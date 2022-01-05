#include <chrono>
#include <fstream>
#include <random>
// #include<function>

#include "GaussNewtonInversion.h"
#include "timer.h"
void func1();
void func2(double);
int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "A thresholding parameter should follow the program name. "
                "Usage:\n";
        cout << "program_name [relative threshold]" << endl;
        return 1;
    }
    double relative_threshold = atof(argv[1]);
    func1();
    func2(relative_threshold);
    return 0;
}
void func1(){
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

    tmesh.out_model_vtk("test_model.vtk");
    tmesh.out_model_netcdf("test_model.nc");

    /******************Foward modelling*****************/
    VectorXd rho;
    tmesh.get_model_parameter_from_mesh(rho, 0);
    Observation ob;

    int nlat = 41;
    int nlon = 41;
    double observation_r = 6428137;
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
    timer.stop();
    cout << "Time for calculating sensitivity: " << timer.getElapsedTimeInSec()
         << " s" << endl;
    timer.start();
    Eigen::VectorXd d_obs = G * rho;
    cout << "Calculation of synthetic data completed" << endl;
    timer.stop();
    cout << "Multiplication time: " << timer.getElapsedTimeInSec() << " s"
         << endl;

    VectorXd test_vec(d_obs.size());
    for (int i = 0; i < test_vec.size(); ++i) {
        test_vec(i) = 1;
    }
    timer.start();
    VectorXd GTv1 = G.transpose() * test_vec;
    timer.stop();
    cout << "Multiplication time of GT*v: " << timer.getElapsedTimeInSec() << "s"
         << endl;

    ofstream ofs("GTv1.txt");
    Mesh m = forward.get_mesh();
    for (int i = 0; i < GTv1.size(); i++) {
        int id = m.get_reordered_id(i);
        ofs << setw(23) << scientific << GTv1(id);
        ofs << endl;
    }         

    // cout << "||d_obs||=" << d_obs.norm() << endl;
    // cout << "0.001*d_obs.norm()=" << 0.001 * d_obs.norm() << endl;

    // GaussNewtonInversion inv(tmesh, ob,
    // Compute_T_rr|Compute_T_rtheta|Compute_T_rphi|Compute_T_thetatheta|Compute_T_thetaphi|Compute_T_phiphi);
    GaussNewtonInversion inv(tmesh, ob, Compute_g_r);

    inv.set_dobs(d_obs);
    inv.output_obs_data("dobs");
}

void func2(double relative_threshold){
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

    tmesh.out_model_vtk("test_model.vtk");
    tmesh.out_model_netcdf("test_model.nc");

        VectorXd rho;
    tmesh.get_model_parameter_from_mesh(rho, 0);
    Observation ob;

    int nlat = 41;
    int nlon = 41;
    double observation_r = 6428137;
    ob.generate_geographic(1, observation_r, observation_r, nlat, lats[0],
                           lats[1], nlon, lons[0], lons[1]);
    ofstream os("sites");
    os << ob;

        //***********************************WAVELET
    Fwd forward_wl(tmesh, ob, Compute_g_r, 4);
    forward_wl.set_use_wavelet(true);
    Timer timer2;
    timer2.start();
    forward_wl.set_compression_threshold(relative_threshold);
    cout << endl;
    cout
        << "=====================Using wavelet transfrom======================="
        << endl;
    cout << "Using wavelet transform\nGenerating synthetic data..." << endl;
    forward_wl.compute_G_wavelet();
    timer2.stop();
    cout << "Time for calculating sensitivity: " << timer2.getElapsedTimeInSec()
         << " s" << endl;

    VectorXd d_obs_wl;
    timer2.start();
    forward_wl.G_vec_mul(rho, d_obs_wl);
    timer2.stop();
    cout << "Multiplication time: " << timer2.getElapsedTimeInSec() << " s"
         << endl;

    VectorXd test_vec(ob.get_n_obs());
    for (int i = 0; i < test_vec.size(); ++i) {
        test_vec(i) = 1;
    }
    timer2.start();
    VectorXd GTv2;
    forward_wl.GT_vec_mul(test_vec, GTv2);
    timer2.stop();
    cout << "Multiplication time of GT*v using wavelet compression: "
         << timer2.getElapsedTimeInSec() << endl;

    ofstream ofs("GTv2.txt");
    Mesh m = forward_wl.get_mesh();
    for (int i = 0; i < GTv2.size(); i++) {
        int id = m.get_reordered_id(i);
        ofs << setw(23) << scientific << GTv2(id);
        ofs << endl;
    }

    GaussNewtonInversion inv2(tmesh, ob, Compute_g_r);
    inv2.set_dobs(d_obs_wl);
    inv2.output_obs_data("dobs_wavelet");
}
