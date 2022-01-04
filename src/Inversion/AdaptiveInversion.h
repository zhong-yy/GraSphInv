#ifndef ADAPTIVE_INVERSION_H
#define ADAPTIVE_INVERSION_H
#include <algorithm>
#include <fstream>

#include "InversionBase.h"

// AdaptiveInversion--Adaptive inversion interface
class AdaptiveInversion : public InversionBase {
   protected:
    double
        refinement_percentage;  // percentage of cells to be refined every time
    int max_refinement_number;  // maximal times of refinement
    int interval_between_refinements;  // interval of iterations beween
                                       // successive refinements
    double min_size_r;
    double min_size_lat;
    double min_size_lon;

   public:
    AdaptiveInversion();
    AdaptiveInversion(const Mesh& mesh_, const Observation& ob_,
                      unsigned long long field_flag_);
    virtual ~AdaptiveInversion();

    void set_max_refinement_number(int n) { max_refinement_number = n; }

    void set_refinement_percentage(double x) { refinement_percentage = x; }

    void set_min_cell_size_in_adaptive_mesh(double min_lat, double min_lon,
                                            double min_r) {
        this->min_size_r = min_r;
        this->min_size_lat = min_lat;
        this->min_size_lon = min_lon;
    }

    void set_interval_between_refinements(int n) {
        interval_between_refinements = n;
    }

    // update sensitivity matrix of gravity forward modelling after mesh is
    // refined
    void expand_G(const map<unsigned int, Cell*>& split_cells);

    // calculate refinement index for each grid cell
    void indicator_calculator(VectorXd& indicator);

    void refine_mesh(double a);

    /**
     * @brief refine inversion mesh and constraint model mesh. The constraint
     * model values in the updated mesh are interpolated from the initial
     * pre-defined interpolator which has been built from gridded data from
     * file.
     *
     * @param  a                  refinement percentage
     * @param  interp             interpolator
     * @param  reference_surface  radius of the reference shpere
     * @param  flag               It must be "crg" or "pet", specifying
     * cross-gradient model or reference density model is interpolated for a
     * predefined "interpolator"
     *
     * see: InversionBase class, InversionBase::create_crg_model_from_data(...)
     *      InversionBase::create_ref_model_from_data(...)
     */
    void refine_mesh(double a, InterpMultilinear<3, double>& interp,
                     double reference_surface = 6378137, string flag = "crg");

    void refine_mesh(double a, InterpMultilinear<3, double>& interp_m0s,
                     InterpMultilinear<3, double>& interp_m0,
                     double reference_surface = 6378137);

    void sort_vec(const VectorXd& vec, VectorXd& sorted_vec, VectorXi& ind) {
        //[0 1 2 3 ... N-1]
        ind = VectorXi::LinSpaced(vec.size(), 0, vec.size() - 1);

        auto rule = [vec](int i, int j) -> bool {
            return vec(i) > vec(j);
        };  //正则表达式，作为sort的谓词
        std::sort(ind.data(), ind.data() + ind.size(), rule);
        // data成员函数返回VectorXd的第一个元素的指针，类似于begin()
        sorted_vec.resize(vec.size());
        for (int i = 0; i < vec.size(); i++) {
            sorted_vec(i) = vec(ind(i));
        }
    }
};

#endif
