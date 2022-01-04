
#ifndef _GUASS_3D_H
#define _GUASS_3D_H
// C++ includes
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

// Local includes
#include "point.h"

// -----------------------------------------------------------------------
// 1D Gauss class definition [-1,+1] segment
class Gauss1D
{
 public:
  Gauss1D (const unsigned int p); // p-point rule for 2p-1 order.
  ~Gauss1D() {}  
  void init (const unsigned int p);  

  unsigned int n_points() { return _points.size();}
  const std::vector<double>&  get_points()  { return _points;}
  const std::vector<double>&  get_weights() { return _weights;} 

  friend std::ostream& operator <<(std::ostream& os, const Gauss1D& Gauss1D);
  
 protected:  
  std::vector<double>            _points;
  std::vector<double>            _weights;  
};

// -----------------------------------------------------------------------
// 2D Gauss class definition [-1,1]x[-1,1] rectangle
class Gauss2D
{
 public:
  Gauss2D (const unsigned int p); // p-point rule for 2p-1 order.
  Gauss2D () {}
  ~Gauss2D() {}  
  void init (const unsigned int p);  

  unsigned int n_points() { return _points.size();}
  const std::vector<Point>&  get_points()  { return _points;}
  const std::vector<double>&   get_weights() { return _weights;} 

  friend std::ostream& operator <<(std::ostream& os, const Gauss2D& Gauss2D);
  
 protected:  
  std::vector<Point>           _points;
  std::vector<double>            _weights;  
};



// -----------------------------------------------------------------------
// 3D Gauss class definition [-1,1]x[-1,1]x[-1,1] cube
class Gauss3D
{
 public:
  Gauss3D (const unsigned int p);
  Gauss3D () {}
  ~Gauss3D() {}  
  void init (const unsigned int p);  

  unsigned int n_points() { return _points.size();}
  const std::vector<Point>&    get_points()  { return _points;}
  const std::vector<double>&   get_weights() { return _weights;} 

  friend std::ostream& operator <<(std::ostream& os, const Gauss3D& Gauss3D);
  
 protected:  
  std::vector<Point>           _points;
  std::vector<double>          _weights;  
};

#endif // _GUASS_3D_H
