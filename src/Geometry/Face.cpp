#include "Face.h"
Face::Face() {
  rc = 0;
  thetac = 0;
  phic = 0;
  for (int i = 0; i < 2; i++) {
    neigh_cells[i] = NULL;
  }
  for (int i = 0; i < 4; i++) {
    child_faces[i] = NULL;
  }
  level = 0;
  isleaf = true;
}
Face::Face(Cell* c1, Cell* c2, unsigned int direction, int level_, bool leaf) {
  this->direction = direction;
  rc = 0;
  thetac = 0;
  phic = 0;
  this->set_neigh_cells(c1, c2, direction);
  this->level = level_;
  isleaf = leaf;
  for (int i = 0; i < 4; i++) {
    child_faces[i] = NULL;
  }
  // if (direction == RADIUS)
  // {
  //     if (c1 != NULL && c2 != NULL)
  //     {
  //         if (!(abs(neigh_cells[0]->_r[1] - rc) < 1e-10))
  //         {
  //             this->display();
  //         }
  //         this->display();
  //         assert(abs(neigh_cells[0]->_r[1] - neigh_cells[1]->_r[0]) < 1e-9);
  //         assert(abs(neigh_cells[0]->_r[1] - rc) < 1e-9);
  //     }
  // }
}
// Face::Face(Cell *c1, Cell *c2, double r, double theta, double phi, int
// level_, bool leaf)
// {
//     rc = r;
//     thetac = theta;
//     phic = phi;
//     this->neigh_cells[0] = c1;
//     this->neigh_cells[1] = c2;
//     this->level = level_;
//     isleaf = leaf;
// }
Face::Face(Cell* c1,
           Cell* c2,
           double r,
           double theta,
           double phi,
           unsigned int direction,
           int level_,
           bool leaf) {
  this->direction = direction;
  rc = r;
  thetac = theta;
  phic = phi;
  this->set_neigh_cells(c1, c2, direction);
  this->level = level_;
  isleaf = leaf;
  if (direction == RADIUS) {
    if (c1 != NULL && c2 != NULL) {
      if (!(abs(neigh_cells[0]->_r[1] - rc) < 1e-10)) {
        display();
      }
      assert(abs(neigh_cells[0]->_r[1] - neigh_cells[1]->_r[0]) < 1e-10);
      assert(abs(neigh_cells[0]->_r[1] - rc) < 1e-10);
    }
  }
  for (int i = 0; i < 4; i++) {
    child_faces[i] = NULL;
  }
}

void Face::set_neigh_cells(Cell* c1, Cell* c2, unsigned int normal_direction) {
  this->neigh_cells[0] = c1;
  this->neigh_cells[1] = c2;
  double r1, theta1, phi1;
  double r2, theta2, phi2;
  if (c1 != NULL && c2 != NULL) {
    c1->get_center(r1, theta1, phi1);
    c2->get_center(r2, theta2, phi2);
    if (normal_direction == NORTH_SOUTH) {
      if (theta1 > theta2) {
        neigh_cells[0] = c2;
        neigh_cells[1] = c1;
        neigh_cells[0]->get_center(r1, theta1, phi1);
        neigh_cells[1]->get_center(r2, theta2, phi2);
        assert(theta1 < theta2);
      }
    } else if (normal_direction == WEST_EAST) {
      if (phi1 > phi2) {
        neigh_cells[0] = c2;
        neigh_cells[1] = c1;

        neigh_cells[0]->get_center(r1, theta1, phi1);
        neigh_cells[1]->get_center(r2, theta2, phi2);
        assert(phi1 < phi2);
      }
    } else if (normal_direction == RADIUS) {
      if (r1 > r2) {
        neigh_cells[0] = c2;
        neigh_cells[1] = c1;

        neigh_cells[0]->get_center(r1, theta1, phi1);
        neigh_cells[1]->get_center(r2, theta2, phi2);
        assert(r1 < r2);
      }
    } else {
      cout << normal_direction << endl;
      this->display();
      cerr << "not NORTH_SOUTH/WEST_EAST/RADIUS" << endl;
      abort();
    }
  }
}

void Face::set_child_faces(Face* f1, Face* f2, Face* f3, Face* f4) {
  this->child_faces[0] = f1;
  this->child_faces[1] = f2;
  this->child_faces[2] = f3;
  this->child_faces[3] = f4;
}
void Face::display() {
  double x, y, z;
  double dx, dy, dz;
  for (int i = 0; i < 2; i++) {
    cout << "Side " << i << ":\t";
    if (neigh_cells[i] == NULL) {
      cout << "NULL" << endl;
    } else {
      cout << "[(" << setprecision(10) << neigh_cells[i]->_r[0] << ", "
           << setprecision(10) << neigh_cells[i]->_r[1] << ")  "
           << "(" << neigh_cells[i]->_theta[0] * 180.0 / GS::PI << ", "
           << neigh_cells[i]->_theta[1] * 180.0 / GS::PI << ")  "
           << "(" << neigh_cells[i]->_phi[0] * 180.0 / GS::PI << ", "
           << neigh_cells[i]->_phi[1] * 180.0 / GS::PI << ")]" << endl;
    }
  }
  cout << "center: (" << rc << ", " << thetac * 180.0 / GS::PI << ", "
       << phic * 180.0 / GS::PI << ")" << endl
       << endl;
}

void Face::set_center(double r, double theta, double phi) {
  this->rc = r;
  this->thetac = theta;
  this->phic = phi;
}