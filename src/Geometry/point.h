
#ifndef _POINT_H
#define _POINT_H

// C++ includes
#include <cmath>
#include <cassert>
#include <iostream>

class Point
{
public:
  Point(const double x = 0., const double y = 0., const double z = 0.);
  Point(const Point &p);
  virtual ~Point();

  friend std::ostream &operator<<(std::ostream &os, const Point &point);

  // Overload operators.
  Point &operator=(const Point &p);
  double operator()(const unsigned int i) const;
  double &operator()(const unsigned int i);
  Point operator+(const Point &v) const;
  Point operator-(const Point &v) const;
  Point operator*(const double a) const;
  double operator*(const Point &v) const;
  bool operator==(const Point &v) const;
  bool operator<(const Point &v) const;
  void operator=(const double a);

  // Math functions
  Point cross(const Point &v) const;
  Point unit() const;
  double size() const;
  void zero();
  double *get_xyz() { return _coords; }
  void set_xyz(double x, double y, double z);

protected:
  double _coords[3];
  const double TOL;
};

//-----------------------------------------------------------------------
inline Point::Point(const double x,
                    const double y,
                    const double z) : TOL(1.0e-10)
{
  _coords[0] = x;
  _coords[1] = y;
  _coords[2] = z;
}

inline Point::Point(const Point &p) : TOL(1.0e-10)
{
  for (unsigned int i = 0; i < 3; i++)
    _coords[i] = p._coords[i];
}

inline Point &Point::operator=(const Point &p)
{
  _coords[0] = p._coords[0];
  _coords[1] = p._coords[1];
  _coords[2] = p._coords[2];

  return (*this);
}

inline Point::~Point()
{
  //no space is allocated by the new operator.
  //so,there is no delete [].
}

inline void Point::set_xyz(double x, double y, double z)
{
  this->_coords[0] = x;
  this->_coords[1] = y;
  this->_coords[2] = z;
}

inline double Point::operator()(const unsigned int i) const
{
  assert(i < 3);
  return _coords[i];
}

inline double &Point::operator()(const unsigned int i)
{
  assert(i < 3);
  return _coords[i];
}

inline Point Point::operator+(const Point &p) const
{
  return Point(_coords[0] + p._coords[0],
               _coords[1] + p._coords[1],
               _coords[2] + p._coords[2]);
}

inline Point Point::operator-(const Point &p) const
{
  return Point(_coords[0] - p._coords[0],
               _coords[1] - p._coords[1],
               _coords[2] - p._coords[2]);
}

inline Point Point::operator*(const double factor) const
{
  return Point(_coords[0] * factor, _coords[1] * factor, _coords[2] * factor);
}

inline double Point::operator*(const Point &p) const
{
  return (_coords[0] * p(0) + _coords[1] * p(1) + _coords[2] * p(2));
}

inline bool Point::operator==(const Point &rhs) const
{
  return ((std::abs(_coords[0] - rhs._coords[0]) +
           std::abs(_coords[1] - rhs._coords[1]) +
           std::abs(_coords[2] - rhs._coords[2])) < 3 * TOL);
}

inline bool Point::operator<(const Point &rhs) const
{
  //First we assume (this)<rhs true
  if (*this == rhs)
    return false;
  if ((*this)(0) < rhs(0))
    return true; //  <
  else if ((*this)(0) > rhs(0))
    return false; //  >
  else if (std::abs((*this)(0) - rhs(0)) < TOL)
  { //vx=rhsx
    if ((*this)(1) < rhs(1))
      return true;
    else if ((*this)(1) > rhs(1))
      return false;
    else if (std::abs((*this)(1) - rhs(1)) < TOL)
    { //vy=rhsy
      if ((*this)(2) < rhs(2))
        return true;
      else if ((*this)(2) > rhs(2))
        return false;
    }
  }
  return false;
}

inline void Point::operator=(const double a)
{
  _coords[0] = a;
  _coords[1] = a;
  _coords[2] = a;
}

inline Point Point::cross(const Point &p) const
{
  return Point(_coords[1] * p._coords[2] - _coords[2] * p._coords[1],
               -_coords[0] * p._coords[2] + _coords[2] * p._coords[0],
               _coords[0] * p._coords[1] - _coords[1] * p._coords[0]);
}

inline Point Point::unit() const
{
  const double length = size();
  return Point(_coords[0] / length,
               _coords[1] / length,
               _coords[2] / length);
}

inline double Point::size() const
{
  double value = std::pow(_coords[0], 2) +
                 std::pow(_coords[1], 2) +
                 std::pow(_coords[2], 2);

  return std::sqrt(value);
}

inline void Point::zero()
{
  _coords[0] = 0.;
  _coords[1] = 0.;
  _coords[2] = 0.;
}

inline std::ostream &operator<<(std::ostream &os, const Point &point)
{
  os << point(0) << "\t" << point(1) << "\t" << point(2);
  return os;
}

#endif // _POINT_H
