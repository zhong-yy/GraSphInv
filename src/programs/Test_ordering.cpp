#include "GaussNewtonInversion.h"
#include "Mesh.h"

int main() {
    double lats[2] = {20, 50};
    double lons[2] = {90, 130};
    double dep[2] = {0, 400000};
    double reference_surface = 6378137;

    Mesh tmesh;
    tmesh.generate_regular_geographic_mesh(lats, 3, lons, 4, dep, 4,
                                           reference_surface);
    tmesh.show_ordering();
    tmesh.set_block_parameter(1, 1, 1, 2, 1, 1, 100);
    VectorXd rho;
    tmesh.get_model_parameter_from_mesh(rho, 0);
    tmesh.out_model_vtk("test_ordering_model.vtk");

    // test ordering after refinement
    int nlat = 3;
    int nlon = 3;
    double observation_r = 6428137;
    Observation ob;
    ob.generate_geographic(1, observation_r, observation_r, nlat, lats[0],
                           lats[1], nlon, lons[0], lons[1]);
    GaussNewtonInversion inv(tmesh, ob, Compute_g_r);
    inv.compute_G();
    inv.set_m(rho);
    inv.set_density_to_mesh();
    inv.refine_mesh(0.05);
    inv.result2vtk("test_ordering_model_refined");
    Mesh& invmesh = inv.get_mesh();
    invmesh.show_ordering();

    return 0;
}
