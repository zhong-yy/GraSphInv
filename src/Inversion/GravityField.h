#ifndef _GRAVITY_FIELD_H
#define _GRAVITY_FIELD_H
// an interface to package "tesseroids-1.2.0"
// http://tesseroids.leouieda.com/en/latest/install.html

// C++ includes
#include <bitset>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
// Tesseroid C library

#include "gauss.h"
#include "gs.h"
#include "tesseroid.h"

// this class includes forward modelling interface to open-source program
// TESSEROIDS, and forward modelling codes of our own.
class GravityField {
   public:
    GravityField(int order_);
    ~GravityField();

    void field_for_a_tesseroid(const Point &x, Tesseroid *tess, double rho,
                               std::vector<double> &field,
                               std::bitset<10> flag = std::bitset<10>(~0ULL));

    void field_for_a_tesseroid_surf(
        const Point &x, Tesseroid *tess, double rho, std::vector<double> &field,
        std::bitset<10> flag = std::bitset<10>(~0ULL));

    // TODO:check
    void field_for_a_tesseroid_adapt(
        const Point &x, Tesseroid *tess, double rho, std::vector<double> &field,
        std::bitset<10> flag = std::bitset<10>(~0ULL));


   private:
    // Gauss Rule used by the Tesseroid library
    int order;
    Gauss3D *unit_cube_gauss;
    Gauss2D *unit_gauss_2D;
    Gauss1D *unit_gauss_1D;
    unsigned int size;

    double kernel_V(const double &l, const double delta[3]) {
        return (1.0 / l);
    }
    double kernel_gr(const double &l, const double delta[3]) {
        double kernel = delta[0] / (l * l * l);
        return kernel;
    };
    double kernel_gtheta(const double &l, const double delta[3]) {
        double kernel = delta[1] / (l * l * l);
        return kernel;
    };
    double kernel_gphi(const double &l, const double delta[3]) {
        double kernel = delta[2] / (l * l * l);
        return kernel;
    }
    double kernel_T_rr(const double &l, const double delta[3]) {
        double l3 = 1.0 / (l * l * l);
        double l2 = 1.0 / (l * l);
        double deltaideltaj = delta[0] * delta[0];
        double kernel = l3 * (3 * deltaideltaj * l2 - 1.0);
        return kernel;
    };
    double kernel_T_rtheta(const double &l, const double delta[3]) {
        double l3 = 1.0 / (l * l * l);
        double l2 = 1.0 / (l * l);
        double deltaideltaj = delta[0] * delta[1];
        double kernel = l3 * (3 * deltaideltaj * l2);
        return kernel;
    }
    double kernel_T_rphi(const double &l, const double delta[3]) {
        double l3 = 1.0 / (l * l * l);
        double l2 = 1.0 / (l * l);
        double deltaideltaj = delta[0] * delta[2];
        double kernel = l3 * (3 * deltaideltaj * l2);
        return kernel;
    }
    double kernel_T_thetatheta(const double &l, const double delta[3]) {
        double l3 = 1.0 / (l * l * l);
        double l2 = 1.0 / (l * l);
        double deltaideltaj = delta[1] * delta[1];
        double kernel = l3 * (3.0 * deltaideltaj * l2 - 1.0);
        return kernel;
    }
    double kernel_T_thetaphi(const double &l, const double delta[3]) {
        double l3 = 1.0 / (l * l * l);
        double l2 = 1.0 / (l * l);
        double deltaideltaj = delta[1] * delta[2];
        double kernel = l3 * (3 * deltaideltaj * l2);
        return kernel;
    }
    double kernel_T_phiphi(const double &l, const double delta[3]) {
        double l3 = 1.0 / (l * l * l);
        double l2 = 1.0 / (l * l);
        double deltaideltaj = delta[2] * delta[2];
        double kernel = l3 * (3.0 * deltaideltaj * l2 - 1.0);
        return kernel;
    }

    double kernel_V_S(const double &r, const double &theta, const double &phi,
                      const double &ro, const double &thetao,
                      const double &phio, int flag) {
        double temp = 0.0, ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);
        if (flag == 1) {
            temp = (r - ro * (sin(theta) * sin(thetao) * cos(phio - phi) +
                              cos(thetao) * cos(theta))) *
                   (r * r * sin(theta)) / ell;
        } else if (flag == 2) {
            temp = ro *
                   (cos(thetao) * sin(theta) -
                    sin(thetao) * cos(theta) * cos(phio - phi)) *
                   sin(theta) * kernel_I_temp(r, ro, ell, cosPsi);
        } else if (flag == 3) {
            temp = ro * sin(thetao) * sin(phi - phio) *
                   kernel_I_temp(r, ro, ell, cosPsi);
        }
        return temp;
    }
    double kernel_gr_S(const double &r, const double &theta, const double &phi,
                       const double &ro, const double &thetao,
                       const double &phio, int flag) {
        Point erprime(sin(thetao) * cos(phio), sin(thetao) * sin(phio),
                      cos(thetao));
        Point temp(0, 0, 0);
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        if (flag == 1) {
            double a = r * r * sin(theta) / ell;
            temp(0) = a * sin(theta) * cos(phi);
            temp(1) = a * sin(theta) * sin(phi);
            temp(2) = a * cos(theta);
        } else if (flag == 2) {
            double a = sin(theta) * kernel_I_temp(r, ro, ell, cosPsi);
            temp(0) = cos(theta) * cos(phi) * a;
            temp(1) = cos(theta) * sin(phi) * a;
            temp(2) = -sin(theta) * a;
        } else if (flag == 3) {
            double a = kernel_I_temp(r, ro, ell, cosPsi);
            temp(0) = -sin(phi) * a;
            temp(1) = cos(phi) * a;
        }
        double result = temp * erprime;
        return result;
    }
    double kernel_gtheta_S(const double &r, const double &theta,
                           const double &phi, const double &ro,
                           const double &thetao, const double &phio, int flag) {
        Point ethetaprime(cos(thetao) * cos(phio), cos(thetao) * sin(phio),
                          -sin(thetao));
        Point temp(0, 0, 0);
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        if (flag == 1) {
            double a = r * r * sin(theta) / ell;
            temp(0) = a * sin(theta) * cos(phi);
            temp(1) = a * sin(theta) * sin(phi);
            temp(2) = a * cos(theta);
        } else if (flag == 2) {
            double a = sin(theta) * kernel_I_temp(r, ro, ell, cosPsi);
            temp(0) = cos(theta) * cos(phi) * a;
            temp(1) = cos(theta) * sin(phi) * a;
            temp(2) = -sin(theta) * a;
        } else if (flag == 3) {
            double a = kernel_I_temp(r, ro, ell, cosPsi);
            temp(0) = -sin(phi) * a;
            temp(1) = cos(phi) * a;
        }
        double result = temp * ethetaprime;
        return result;
    }
    double kernel_gphi_S(const double &r, const double &theta,
                         const double &phi, const double &ro,
                         const double &thetao, const double &phio, int flag) {
        Point ephiprime(-sin(phio), cos(phio), 0);
        Point temp(0, 0, 0);
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        if (flag == 1) {
            double a = r * r * sin(theta) / ell;
            temp(0) = a * sin(theta) * cos(phi);
            temp(1) = a * sin(theta) * sin(phi);
            temp(2) = a * cos(theta);
        } else if (flag == 2) {
            double a = sin(theta) * kernel_I_temp(r, ro, ell, cosPsi);
            temp(0) = cos(theta) * cos(phi) * a;
            temp(1) = cos(theta) * sin(phi) * a;
            temp(2) = -sin(theta) * a;
        } else if (flag == 3) {
            double a = kernel_I_temp(r, ro, ell, cosPsi);
            temp(0) = -sin(phi) * a;
            temp(1) = cos(phi) * a;
        }
        double result = temp * ephiprime;
        return result;
    }
    double kernel_Trr_S(const double &r, const double &theta, const double &phi,
                        const double &ro, const double &thetao,
                        const double &phio, int flag) {
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        Point e1(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
        Point e2(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
        Point e3(-sin(phi), cos(phi));

        Point e1prime(cos(phio) * sin(thetao), sin(phio) * sin(thetao),
                      cos(thetao));
        Point e2prime(cos(thetao) * cos(phio), cos(thetao) * sin(phio),
                      -sin(thetao));
        Point e3prime(-sin(phio), cos(phio));
        double result = 0;
        double A1 = e1 * e1prime;
        double B1 = e2 * e1prime;
        double C1 = e3 * e1prime;
        if (flag == 1) {
            result =
                (ro - r * A1) * A1 * r * r * sin(theta) / (ell * ell * ell);
        } else if (flag == 2) {
            result = (ro * kernel_J_temp(r, ro, ell, cosPsi) -
                      A1 * kernel_K_temp(r, ro, ell, cosPsi)) *
                     B1 * sin(theta);
        } else if (flag == 3) {
            result = (ro * kernel_J_temp(r, ro, ell, cosPsi) -
                      A1 * kernel_K_temp(r, ro, ell, cosPsi)) *
                     C1;
        }
        return result;
    }
    double kernel_Trtheta_S(const double &r, const double &theta,
                            const double &phi, const double &ro,
                            const double &thetao, const double &phio,
                            int flag) {
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        Point e1(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
        Point e2(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
        Point e3(-sin(phi), cos(phi));

        Point e1prime(cos(phio) * sin(thetao), sin(phio) * sin(thetao),
                      cos(thetao));
        Point e2prime(cos(thetao) * cos(phio), cos(thetao) * sin(phio),
                      -sin(thetao));
        Point e3prime(-sin(phio), cos(phio));
        double result = 0;
        double A1 = e1 * e1prime;
        double A2 = e1 * e2prime;
        double B2 = e2 * e2prime;
        double C2 = e3 * e2prime;
        if (flag == 1) {
            result =
                (ro - r * A1) * A2 * r * r * sin(theta) / (ell * ell * ell);
        } else if (flag == 2) {
            result = (ro * kernel_J_temp(r, ro, ell, cosPsi) -
                      A1 * kernel_K_temp(r, ro, ell, cosPsi)) *
                     B2 * sin(theta);
        } else if (flag == 3) {
            result = (ro * kernel_J_temp(r, ro, ell, cosPsi) -
                      A1 * kernel_K_temp(r, ro, ell, cosPsi)) *
                     C2;
        }
        return result;
    }
    double kernel_Trphi_S(const double &r, const double &theta,
                          const double &phi, const double &ro,
                          const double &thetao, const double &phio, int flag) {
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        Point e1(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
        Point e2(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
        Point e3(-sin(phi), cos(phi));

        Point e1prime(cos(phio) * sin(thetao), sin(phio) * sin(thetao),
                      cos(thetao));
        Point e2prime(cos(thetao) * cos(phio), cos(thetao) * sin(phio),
                      -sin(thetao));
        Point e3prime(-sin(phio), cos(phio));
        double result = 0;
        double A1 = e1 * e1prime;
        double A3 = e1 * e3prime;
        double B3 = e2 * e3prime;
        double C3 = e3 * e3prime;
        if (flag == 1) {
            result =
                (ro - r * A1) * A3 * r * r * sin(theta) / (ell * ell * ell);
        } else if (flag == 2) {
            result = (ro * kernel_J_temp(r, ro, ell, cosPsi) -
                      A1 * kernel_K_temp(r, ro, ell, cosPsi)) *
                     B3 * sin(theta);
        } else if (flag == 3) {
            result = (ro * kernel_J_temp(r, ro, ell, cosPsi) -
                      A1 * kernel_K_temp(r, ro, ell, cosPsi)) *
                     C3;
        }
        return result;
    }
    double kernel_Tthetatheta_S(const double &r, const double &theta,
                                const double &phi, const double &ro,
                                const double &thetao, const double &phio,
                                int flag) {
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        Point e1(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
        Point e2(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
        Point e3(-sin(phi), cos(phi));

        Point e1prime(cos(phio) * sin(thetao), sin(phio) * sin(thetao),
                      cos(thetao));
        Point e2prime(cos(thetao) * cos(phio), cos(thetao) * sin(phio),
                      -sin(thetao));
        Point e3prime(-sin(phio), cos(phio));
        double result = 0;
        double A2 = e1 * e2prime;
        double B2 = e2 * e2prime;
        double C2 = e3 * e2prime;
        if (flag == 1) {
            result = -r * A2 * A2 * r * r * sin(theta) / (ell * ell * ell);
        } else if (flag == 2) {
            result = -A2 * kernel_K_temp(r, ro, ell, cosPsi) * B2 * sin(theta);
        } else if (flag == 3) {
            result = -A2 * kernel_K_temp(r, ro, ell, cosPsi) * C2;
        }
        return result;
    }
    double kernel_Tthetaphi_S(const double &r, const double &theta,
                              const double &phi, const double &ro,
                              const double &thetao, const double &phio,
                              int flag) {
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        Point e1(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
        Point e2(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
        Point e3(-sin(phi), cos(phi));

        Point e1prime(cos(phio) * sin(thetao), sin(phio) * sin(thetao),
                      cos(thetao));
        Point e2prime(cos(thetao) * cos(phio), cos(thetao) * sin(phio),
                      -sin(thetao));
        Point e3prime(-sin(phio), cos(phio));
        double result = 0;
        double A2 = e1 * e2prime;
        double A3 = e1 * e3prime;
        double B3 = e2 * e3prime;
        double C3 = e3 * e3prime;
        if (flag == 1) {
            result = -r * A2 * A3 * r * r * sin(theta) / (ell * ell * ell);
        } else if (flag == 2) {
            result = -A2 * kernel_K_temp(r, ro, ell, cosPsi) * B3 * sin(theta);
        } else if (flag == 3) {
            result = -A2 * kernel_K_temp(r, ro, ell, cosPsi) * C3;
        }
        return result;
    }
    double kernel_Tphiphi_S(const double &r, const double &theta,
                            const double &phi, const double &ro,
                            const double &thetao, const double &phio,
                            int flag) {
        double ell = 0.0;
        double cosPsi = 0.0;
        cosPsi = cos(theta) * cos(thetao) +
                 sin(theta) * sin(thetao) * cos(phi - phio);
        ell = sqrt(r * r + ro * ro - 2.0 * cosPsi * r * ro);

        Point e1(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
        Point e2(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
        Point e3(-sin(phi), cos(phi));

        Point e1prime(cos(phio) * sin(thetao), sin(phio) * sin(thetao),
                      cos(thetao));
        Point e2prime(cos(thetao) * cos(phio), cos(thetao) * sin(phio),
                      -sin(thetao));
        Point e3prime(-sin(phio), cos(phio));
        double result = 0;
        double A3 = e1 * e3prime;
        double B3 = e2 * e3prime;
        double C3 = e3 * e3prime;
        if (flag == 1) {
            result = -r * A3 * A3 * r * r * sin(theta) / (ell * ell * ell);
        } else if (flag == 2) {
            result = -A3 * kernel_K_temp(r, ro, ell, cosPsi) * B3 * sin(theta);
        } else if (flag == 3) {
            result = -A3 * kernel_K_temp(r, ro, ell, cosPsi) * C3;
        }
        return result;
    }
    double kernel_I_temp(const double &r, const double &ro, const double &ell,
                         const double &cosPsi) {
        double I;
        if (fabs(cosPsi - 1) < TOL) {
            if (ro < r) {
                I = r - ro + ro * log(r - ro);
            } else {
                I = ro - r - ro * log(ro - r);
            }
        } else if (fabs(cosPsi + 1) < TOL) {
            I = r + ro - ro * log(r + ro);
        } else {
            I = ell + ro * cosPsi * log(ell + r - ro * cosPsi);
        }
        return I;
    }
    double kernel_J_temp(const double &r, const double &ro, const double &ell,
                         const double &cosPsi) {
        double J;
        if (fabs(cosPsi - 1) < TOL) {
            if (ro < r) {
                J = (ro - 2.0 * r) / (2.0 * (r - ro) * (r - ro));
            } else {
                J = (2.0 * r - ro) / (2.0 * (r - ro) * (r - ro));
            }
        } else if (fabs(cosPsi + 1) < TOL) {
            J = (-ro - 2.0 * r) / (2.0 * (r + ro) * (r + ro));
        } else {
            J = (r * cosPsi - ro) / (ro * ell * (1 - cosPsi * cosPsi));
        }
        return J;
    }
    double kernel_K_temp(const double &r, const double &ro, const double &ell,
                         const double &cosPsi) {
        double K;
        if (fabs(cosPsi - 1) < TOL) {
            if (ro < r) {
                K = log(r - ro) - 2.0 * ro / (r - ro) -
                    ro * ro / (2.0 * (r - ro) * (r - ro));
            } else {
                K = -(log(ro - r) - 2.0 * ro / (r - ro) -
                      ro * ro / (2.0 * (r - ro) * (r - ro)));
            }
        } else if (fabs(cosPsi + 1) < TOL) {
            K = log(r + ro) + 2.0 * ro / (r + ro) -
                ro * ro / (2.0 * (r + ro) * (r + ro));
        } else {
            K = ((2.0 * cosPsi * cosPsi - 1) * r - ro * cosPsi) /
                    (ell * (1.0 - cosPsi * cosPsi)) +
                log(ell + r - ro * cosPsi);
        }
        return K;
    }
};

#endif  // _GRAVITY_FIELD_H
