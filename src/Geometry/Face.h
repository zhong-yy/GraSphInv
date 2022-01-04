#ifndef _FACE
#define _FACE
#include "Cell.h"
#include<iomanip>
//Face tree
class Cell;
class Face
{
public:
  Face();
  Face(Cell *c1, Cell *c2, unsigned int direction, int level = 0, bool leaf = true);
  // Face(Cell *c1, Cell *c2, double r, double theta, double phi, int level = 0, bool leaf = true);
  Face(Cell *c1, Cell *c2, double r, double theta, double phi, unsigned int direction, int level = 0, bool leaf = true);
  void set_neigh_cells(Cell *c1, Cell *c2, unsigned int normal_direction);
  void set_child_faces(Face *f1, Face *f2, Face *f3, Face *f4);

  void display();

  void set_center(double r, double theta, double phi);
  // protected:
  Cell *neigh_cells[2];
  Face *child_faces[4];

  int level;
  bool isleaf;

  double rc, thetac, phic;
  unsigned int direction;
};

#endif