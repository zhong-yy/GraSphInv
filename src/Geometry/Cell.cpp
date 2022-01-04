#include "Cell.h"
Cell::Cell(double r0, double theta0, double phi0, double r1, double theta1,
           double phi1, int level_, int num_para, bool isleaf_)
    : Tesseroid(r0, theta0, phi0, r1, theta1, phi1) {
    this->id = -1;
    this->parameters.resize(num_para);
    for (int i = 0; i < num_para; i++) {
        this->parameters[i] = 0;
    }

    this->level = level_;
    this->isleaf = isleaf_;

    for (int i = 0; i < 8; i++) {
        child_cells[i] = NULL;
    }
    for (int i = 0; i < 2; i++) {
        external_faces_theta[i] = NULL;
        external_faces_phi[i] = NULL;
        external_faces_r[i] = NULL;
    }
    for (int i = 0; i < 4; i++) {
        internal_faces_theta[i] = NULL;
        internal_faces_phi[i] = NULL;
        internal_faces_r[i] = NULL;
    }
    this->ordering_forward = true;
}

void Cell::set_internal_faces(Face *f[4], unsigned int normal_dirction) {
    if (normal_dirction == NORTH_SOUTH) {
        assert(abs(f[0]->thetac - f[1]->thetac) < TOL &&
               abs(f[0]->thetac - f[2]->thetac) < TOL &&
               abs(f[0]->thetac - f[3]->thetac) < TOL);
        for (int i = 0; i < 3; i++) {
            for (int j = 1; j < 4 - i; j++) {
                if (f[j - 1]->phic > f[j]->phic) {
                    Face *temp = f[j - 1];
                    f[j - 1] = f[j];
                    f[j] = temp;
                }
            }
        }
        if (f[0]->rc > f[1]->rc) {
            Face *temp = f[0];
            f[0] = f[1];
            f[1] = temp;
        }
        if (f[2]->rc > f[3]->rc) {
            Face *temp = f[2];
            f[2] = f[3];
            f[3] = temp;
        }
        for (int i = 0; i < 4; i++) {
            this->internal_faces_theta[i] = f[i];
        }
    } else if (normal_dirction == WEST_EAST) {
        assert(abs(f[0]->phic - f[1]->phic) < TOL &&
               abs(f[0]->phic - f[2]->phic) < TOL &&
               abs(f[0]->phic - f[3]->phic) < TOL);
        for (int i = 0; i < 3; i++) {
            for (int j = 1; j < 4 - i; j++) {
                if (f[j - 1]->thetac > f[j]->thetac) {
                    Face *temp = f[j - 1];
                    f[j - 1] = f[j];
                    f[j] = temp;
                }
            }
        }
        if (f[0]->rc > f[1]->rc) {
            Face *temp = f[0];
            f[0] = f[1];
            f[1] = temp;
        }
        if (f[2]->rc > f[3]->rc) {
            Face *temp = f[2];
            f[2] = f[3];
            f[3] = temp;
        }
        for (int i = 0; i < 4; i++) {
            this->internal_faces_phi[i] = f[i];
        }
    } else if (normal_dirction == RADIUS) {
        assert(abs(f[0]->rc - f[1]->rc) < TOL &&
               abs(f[0]->rc - f[2]->rc) < TOL &&
               abs(f[0]->rc - f[3]->rc) < TOL);
        for (int i = 0; i < 3; i++) {
            for (int j = 1; j < 4 - i; j++) {
                if (f[j - 1]->thetac > f[j]->thetac) {
                    Face *temp = f[j - 1];
                    f[j - 1] = f[j];
                    f[j] = temp;
                }
            }
        }
        if (f[0]->phic > f[1]->phic) {
            Face *temp = f[0];
            f[0] = f[1];
            f[1] = temp;
        }
        if (f[2]->phic > f[3]->phic) {
            Face *temp = f[2];
            f[2] = f[3];
            f[3] = temp;
        }
        for (int i = 0; i < 4; i++) {
            this->internal_faces_r[i] = f[i];
        }
    } else {
        cerr << "Not theta/phi/r" << endl;
    }
}
void Cell::set_external_faces(Face *f1, Face *f2,
                              unsigned int normal_dirction) {
    switch (normal_dirction) {
        case NORTH_SOUTH:
            external_faces_theta[0] = f1;
            external_faces_theta[1] = f2;
            if (f1->thetac > f2->thetac) {
                // cout << "x1>x2, f1: (" << f1->thetac << ", " << f1->phic <<
                // ", " << f1->rc << ")  f2: (" << f2->thetac << ", " <<
                // f2->phic << ", " << f2->rc << ")" << endl;
                external_faces_theta[0] = f2;
                external_faces_theta[1] = f1;
            }
            break;
        case WEST_EAST:
            external_faces_phi[0] = f1;
            external_faces_phi[1] = f2;
            if (f1->phic > f2->phic) {
                external_faces_phi[0] = f2;
                external_faces_phi[1] = f1;
            }
            break;
        case RADIUS:
            external_faces_r[0] = f1;
            external_faces_r[1] = f2;
            if (f1->rc > f2->rc) {
                external_faces_r[0] = f2;
                external_faces_r[1] = f1;
            }
            break;
        default:
            cerr << "Not theta/phi/r" << endl;
    }
}

// void Cell::set_external_faces(int i, Face *face, char normal_dirction)
// {
//     assert(i >= 0 && i <= 1);
//     switch (normal_dirction)
//     {
//     case 'theta':
//         external_faces_theta[i] = face;
//         break;
//     case 'phi':
//         external_faces_phi[i] = face;
//         break;
//     case 'r':
//         external_faces_r[i] = face;
//         break;
//     }
// }
