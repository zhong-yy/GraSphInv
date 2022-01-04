
#ifndef _GS_H
#define _GS_H

#include <complex>
#include <map>
#include <vector>

#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#endif
#include <Eigen/Dense>   // Linear Algebra Lib.
#include <Eigen/Sparse>  // Sparse Lib
#include <Eigen/StdVector>

namespace GS {
static const double PI = 3.1415926535897932384626433832795;
static const unsigned int INVALID_UNIT = static_cast<unsigned int>(-1);
static const double TOL = 1.000e-12;
static const double G0 = 6.673e-11;  // m^3*kg^{-1}*s^{-1}
// static const std::complex<double> I=std::complex<double>(0.0,1.0);
static const double Mean_Earth_Radius = 6378137.0;  // m, 6378.137km
static const double MER = 6378137.0;                // m
static const double SI2Eotvos = 1000000000.0;  // frac{1}{s^2} = 10^9\ Eotvos
static const double SI2mGal = 100000.0;        // frac{m}{s^2} = 10^5\ mGal

typedef std::complex<double> Dcomplex;
typedef std::map<unsigned int, std::vector<unsigned int> > Field;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SMatrix;
typedef Eigen::VectorXd Vector;

typedef Eigen::SparseVector<double> SparseVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
typedef Eigen::Matrix<double, 4, 1> Vector4D;
typedef Eigen::Matrix<double, 4, 4> Matrix4D;
typedef Eigen::Matrix<double, 3, 3> Matrix3D;
typedef Eigen::Matrix<double, 4, 3> Matrix43D;
typedef Eigen::Matrix<double, 1, 3> Matrix13D;
typedef Eigen::Matrix<double, 3, 1> Matrix31D;
// typedef Eigen::Matrix<Dcomplex, 3,   1>    Matrix31C;
// typedef Eigen::Matrix<Dcomplex, 3,   3>    Matrix3C;
typedef Eigen::Matrix<double, 1, 4> Matrix14D;
typedef Eigen::VectorXd DenseVector;
typedef Eigen::Vector3d Vector3D;

enum ComputationOptions {
  Compute_V = 0x0001,
  Compute_g_r = 0x0002,
  Compute_g_theta = 0x0004,
  Compute_g_phi = 0x0008,
  Compute_T_rr = 0x0010,
  Compute_T_rtheta = 0x0020,
  Compute_T_rphi = 0x0040,
  Compute_T_thetatheta = 0x0080,
  Compute_T_thetaphi = 0x0100,
  Compute_T_phiphi = 0x0200
};
enum DIRECTION {
  RADIUS,
  NORTH_SOUTH,
  WEST_EAST,
};
}  // namespace GS

using namespace GS;

//---------------------------------------------
#endif  // #define _GS_H
