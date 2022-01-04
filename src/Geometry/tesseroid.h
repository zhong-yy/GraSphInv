/*
  A Tesseroid is composed of [r1,r2]x[theta1,theta2]x[phi1,phi2]
  ri is the radius [0,infinite]
  thetai, polar angle [0,pi], or colatitude, unit is radian
          0       :  90 N
          pi/2.0  :  equator
          pi      :  90 S
  phii,   azimuthal angle, [0,2pi]
          0       :  180 W
          pi      :  0
          2pi     :  180 E
  r2>r1, theta2>theta1, phi2>phi1
  see https://en.wikipedia.org/wiki/Spherical_coordinate_system
*/

#ifndef _Tesseroid_H
#define _Tesseroid_H

#include "gs.h"
#include "point.h"

class Tesseroid {
   public:
    Tesseroid();
    Tesseroid(unsigned int id, double r[2], double theta[2], double phi[2],
              double RHO);
    Tesseroid(double r0, double theta0, double phi0, double r1, double theta1,
              double phi1);
    ~Tesseroid() {}

    // friend class Tgrid;
    friend std::ostream &operator<<(std::ostream &os, const Tesseroid &T) {
        os << T._id << "\t" << T._r[1] << "\t" << T._r[0] << "\t"
           << T._theta[1] * 180 / GS::PI << "\t" << T._theta[0] * 180 / GS::PI
           << "\t" << T._phi[1] * 180 / GS::PI << "\t"
           << T._phi[0] * 180 / GS::PI << "\t" << T._density << "\n";
        return os;
    }

    unsigned int get_id() const { return _id; }
    // int get_marker() { return _marker; }
    double get_density() const { return _density; }
    void set_density(const double &density) { this->_density = density; }
    double get_size();
    void get_size(double &dr, double &dtheta, double &dphi);
    Point get_gpoint();
    // void set_neighbors(Tesseroid *n, Tesseroid *s, Tesseroid *w, Tesseroid
    // *e, Tesseroid *up, Tesseroid *down);

    double get_volumn() const {
        return ((_r[1] * _r[1] * _r[1] - _r[0] * _r[0] * _r[0]) / 3.0 *
                (cos(_theta[0]) - cos(_theta[1])) * (_phi[1] - _phi[0]));
    }

    void get_center(double &rc, double &thetac, double &phic);

    unsigned int _id;
    double _density;
    double _r[2];
    double _theta[2];  // radians
    double _phi[2];    // radians
    // int _marker;

    // Tesseroid *n_neighbor;    //north neighbour
    // Tesseroid *s_neighbor;    //south
    // Tesseroid *w_neighbor;    //west
    // Tesseroid *e_neighbor;    //east
    // Tesseroid *up_neighbor;   //upside
    // Tesseroid *down_neighbor; //downside
};

inline Tesseroid::Tesseroid() {
    _id = -1;
    for (unsigned int i = 0; i < 2; i++) _r[i] = 0;
    for (unsigned int i = 0; i < 2; i++) _theta[i] = 0;
    for (unsigned int i = 0; i < 2; i++) _phi[i] = 0;
    _density = 0.;
    // _marker = -99;
    // n_neighbor = NULL;
    // s_neighbor = NULL;
    // w_neighbor = NULL;
    // e_neighbor = NULL;
    // up_neighbor = NULL;
    // down_neighbor = NULL;
}

inline Tesseroid::Tesseroid(unsigned int id, double r[2], double theta[2],
                            double phi[2], double rho) {
    _id = id;
    for (unsigned int i = 0; i < 2; i++) _r[i] = r[i];
    for (unsigned int i = 0; i < 2; i++) _theta[i] = theta[i];
    for (unsigned int i = 0; i < 2; i++) _phi[i] = phi[i];
    _density = rho;
    // _marker = marker;

    // n_neighbor = NULL;
    // s_neighbor = NULL;
    // w_neighbor = NULL;
    // e_neighbor = NULL;
    // up_neighbor = NULL;
    // down_neighbor = NULL;
    // check
    assert(_r[1] > _r[0] || std::abs(_r[1] - _r[0]) < TOL);
    assert(_theta[1] > _theta[0]);
    assert(_phi[1] > _phi[0]);
}

inline Tesseroid::Tesseroid(double r0, double theta0, double phi0, double r1,
                            double theta1, double phi1) {
    this->_r[0] = std::min(r0, r1);
    this->_theta[0] = std::min(theta0, theta1);
    this->_phi[0] = std::min(phi0, phi1);

    this->_r[1] = std::max(r0, r1);
    this->_theta[1] = std::max(theta0, theta1);
    this->_phi[1] = std::max(phi0, phi1);
    assert(_r[1] > _r[0] || std::abs(_r[1] - _r[0]) < TOL);
    assert(_theta[1] > _theta[0]);
    assert(_phi[1] > _phi[0]);
    _density = 0;
}

inline void Tesseroid::get_size(double &dr, double &dtheta, double &dphi) {
    dr = this->_r[1] - this->_r[0];
    dtheta = this->_theta[1] - this->_theta[0];
    dphi = this->_phi[1] - this->_phi[0];
}
inline double Tesseroid::get_size() {
    // int_{r1,r2} int_{theta2,theta1} int_{phi1,phi2} J drdthetadphi
    // J = r*sin(theta)
    double vr = (std::pow(_r[1], 3.) - std::pow(_r[0], 3.)) / 3.;
    double vtheta = -1.0 * std::cos(_theta[1]) + std::cos(_theta[0]);
    double vphi = _phi[1] - _phi[0];
    double v = vr * vtheta * vphi;
    assert(v > 0.);
    return v;
}

inline Point Tesseroid::get_gpoint() {
    return Point((_r[1] + _r[0]) * 0.5, (_theta[1] + _theta[0]) * 0.5,
                 (_phi[1] + _phi[0]) * 0.5);
}

inline void Tesseroid::get_center(double &rc, double &thetac, double &phic) {
    rc = 0.5 * (_r[1] + _r[0]);
    thetac = 0.5 * (_theta[1] + _theta[0]);
    phic = 0.5 * (_phi[1] + _phi[0]);
}

// inline void Tesseroid::set_neighbors(Tesseroid *n, Tesseroid *s, Tesseroid
// *w, Tesseroid *e, Tesseroid *up, Tesseroid *down)
// {
//   this->n_neighbor = n;
//   this->s_neighbor = s;
//   this->w_neighbor = w;
//   this->e_neighbor = e;
//   this->up_neighbor = up;
//   this->down_neighbor = down;
//   return;
// }

#endif  // _Tesseroid_h
