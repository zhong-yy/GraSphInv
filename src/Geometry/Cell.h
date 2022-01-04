#ifndef _CELL
#define _CELL
#include <Eigen/Dense>   // Linear Algebra Lib.
#include <Eigen/Sparse>  // Sparse Lib
#include <Eigen/StdVector>
#include <algorithm>
#include <cassert>
#include <string>

#include "Face.h"
#include "tesseroid.h"
using namespace std;
using namespace Eigen;

// Cell tree
class Face;
class Cell : public Tesseroid {
   public:
    // ~Cell();
    Cell(double r0, double theta0, double phi0, double r1, double theta1,
         double phi1, int level_ = 0, int num_para = 1, bool isleaf_ = true);
    bool get_ordering_forward() { return this->ordering_forward; }
    void set_ordering_forward(bool forward) {
        this->ordering_forward = forward;
    }
    void set_id(int _id) { this->id = _id; }
    int get_id() { return id; }
    void set_parameter(double new_value, int i = 0) {
        // cout<<"xxx"<<this->parameters.size()<<endl;
        this->parameters[i] = new_value;
        if (!isleaf) {
            for (int j = 0; j < 8; j++) {
                child_cells[j]->set_parameter(new_value, i);
            }
        }
    }
    double get_parameter(int i = 0) { return this->parameters[i]; }
    void set_external_faces(Face* f1, Face* f2, unsigned int normal_dirction);
    // void set_external_faces(int i, Face *face, char direction);

    void set_internal_faces(Face* f[4], unsigned int normal_dirction);

    int get_level() { return level; }

    //   protected:
    bool isleaf;
    bool ordering_forward;
    int level;  // first level is 0
    int id;
    // Rect_prism *cell;
    Cell* child_cells[8];

    Face* external_faces_r[2];
    Face* external_faces_theta[2];
    Face* external_faces_phi[2];

    Face* internal_faces_r[4];
    Face* internal_faces_theta[4];
    Face* internal_faces_phi[4];

    vector<double> parameters;
    // double density;
};
#endif
