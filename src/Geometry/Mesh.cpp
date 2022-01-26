#include "Mesh.h"
Mesh::Mesh()
    : leaf_cells(0),
      leaf_faces(0),
      ro_leaf_cells(0),
      cells(0),
      faces(0),
      num_cell(0),
      num_face(0),
      new_to_old_index(0),
      old_to_new_index(0) {
    nr = 0;
    ntheta = 0;
    nphi = 0;
    num_leaf_cells = 0;
    num_leaf_faces = 0;
    n_parameters = 1;
    r_lim[0] = 0;
    r_lim[1] = 0;
    theta_lim[0] = 0;
    theta_lim[1] = 0;
    phi_lim[0] = 0;
    phi_lim[1] = 0;
    // cout<<"1"<<endl;
}
Mesh::~Mesh() { this->clear_all(); }
Mesh::Mesh(const Mesh& source_mesh) {
    this->clear_all();  // clear

    // copy non-pointer type data
    this->nr = source_mesh.nr;
    this->ntheta = source_mesh.ntheta;
    this->nphi = source_mesh.nphi;
    for (int i = 0; i < 2; i++) {
        this->r_lim[i] = source_mesh.r_lim[i];
        this->theta_lim[i] = source_mesh.theta_lim[i];
        this->phi_lim[i] = source_mesh.phi_lim[i];
    }
    this->num_face = source_mesh.num_face;
    this->num_cell = source_mesh.num_cell;
    this->r_points = source_mesh.r_points;
    this->num_leaf_cells = source_mesh.num_leaf_cells;
    this->num_leaf_faces = source_mesh.num_leaf_faces;
    this->n_parameters = source_mesh.n_parameters;
    this->new_to_old_index = source_mesh.new_to_old_index;
    this->old_to_new_index = source_mesh.old_to_new_index;

    // cout<<"A"<<endl;
    // deep copy pointers
    unordered_map<Cell*, Cell*> cc_mp;  // key: old mesh, value: new mesh
    unordered_map<Face*, Face*> ff_mp;

    this->cells.resize(source_mesh.cells.size());
    for (int i = 0; i < source_mesh.cells.size(); i++) {
        this->cells[i].resize(source_mesh.cells[i].size());
        for (int j = 0; j < source_mesh.cells[i].size(); j++) {
            Cell* c = source_mesh.cells[i][j];
            this->cells[i][j] = new Cell(
                c->_r[0], c->_theta[0], c->_phi[0], c->_r[1], c->_theta[1],
                c->_phi[1], c->level, c->parameters.size(), c->isleaf);
            this->cells[i][j]->set_id(c->id);
            this->cells[i][j]->_density = c->_density;
            this->cells[i][j]->set_ordering_forward(c->get_ordering_forward());
            cc_mp[c] = this->cells[i][j];
        }
    }

    // cout<<"A"<<endl;

    this->faces.resize(source_mesh.faces.size());
    for (int i = 0; i < source_mesh.faces.size(); i++) {
        this->faces[i].resize(source_mesh.faces[i].size());
        for (int j = 0; j < source_mesh.faces[i].size(); j++) {
            Face* f = source_mesh.faces[i][j];
            this->faces[i][j] = new Face(NULL, NULL, f->rc, f->thetac, f->phic,
                                         f->direction, f->level, f->isleaf);
            ff_mp[f] = this->faces[i][j];
        }
    }

    // set associated pointers
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < cells[i].size(); j++) {
            Cell* c = source_mesh.cells[i][j];
            for (int k = 0; k < 8; k++) {
                if (c->child_cells[k] != NULL) {
                    cc_mp[c]->child_cells[k] = cc_mp[c->child_cells[k]];
                } else {
                    cc_mp[c]->child_cells[k] = NULL;
                }
            }
            // cout<<"X"<<endl;
            for (int k = 0; k < 2; k++) {
                if (c->external_faces_r[k] != NULL) {
                    cc_mp[c]->external_faces_r[k] =
                        ff_mp[c->external_faces_r[k]];
                } else {
                    cc_mp[c]->external_faces_r[k] = NULL;
                }

                if (c->external_faces_theta[k] != NULL) {
                    cc_mp[c]->external_faces_theta[k] =
                        ff_mp[c->external_faces_theta[k]];
                } else {
                    cc_mp[c]->external_faces_theta[k] = NULL;
                }

                if (c->external_faces_phi[k] != NULL) {
                    cc_mp[c]->external_faces_phi[k] =
                        ff_mp[c->external_faces_phi[k]];
                } else {
                    cc_mp[c]->external_faces_phi[k] = NULL;
                }
            }
            for (int k = 0; k < 4; k++) {
                if (c->internal_faces_r[k] != NULL) {
                    cc_mp[c]->internal_faces_r[k] =
                        ff_mp[c->internal_faces_r[k]];
                } else {
                    cc_mp[c]->internal_faces_r[k] = NULL;
                }

                if (c->internal_faces_theta[k] != NULL) {
                    cc_mp[c]->internal_faces_theta[k] =
                        ff_mp[c->internal_faces_theta[k]];
                } else {
                    cc_mp[c]->internal_faces_theta[k] = NULL;
                }

                if (c->internal_faces_phi[k] != NULL) {
                    cc_mp[c]->internal_faces_phi[k] =
                        ff_mp[c->internal_faces_phi[k]];
                } else {
                    cc_mp[c]->internal_faces_phi[k] = NULL;
                }
            }
        }
    }

    for (int i = 0; i < faces.size(); i++) {
        for (int j = 0; j < faces[i].size(); j++) {
            Face* f = source_mesh.faces[i][j];
            for (int k = 0; k < 4; k++) {
                if (f->child_faces[k] != NULL) {
                    ff_mp[f]->child_faces[k] = ff_mp[f->child_faces[k]];
                } else {
                    ff_mp[f]->child_faces[k] = NULL;
                }
            }
            for (int k = 0; k < 2; k++) {
                if (f->neigh_cells[k] != NULL) {
                    ff_mp[f]->neigh_cells[k] = cc_mp[f->neigh_cells[k]];
                } else {
                    ff_mp[f]->neigh_cells[k] = NULL;
                }
            }
        }
    }

    leaf_cells.resize(source_mesh.leaf_cells.size());
    ro_leaf_cells.resize(source_mesh.ro_leaf_cells.size());
    leaf_faces.resize(source_mesh.leaf_faces.size());
    for (int i = 0; i < leaf_cells.size(); i++) {
        Cell* c = source_mesh.leaf_cells[i];
        this->leaf_cells[i] = cc_mp[c];
    }
    for (int i = 0; i < ro_leaf_cells.size(); i++) {
        Cell* c = source_mesh.ro_leaf_cells[i];
        this->ro_leaf_cells[i] = cc_mp[c];
    }
    for (int i = 0; i < leaf_faces.size(); i++) {
        Face* f = source_mesh.leaf_faces[i];
        this->leaf_faces[i] = ff_mp[f];
    }
}

Mesh& Mesh::operator=(const Mesh& source_mesh) {
    this->clear_all();  // clear

    // copy non-pointer type data
    this->nr = source_mesh.nr;
    this->ntheta = source_mesh.ntheta;
    this->nphi = source_mesh.nphi;
    for (int i = 0; i < 2; i++) {
        this->r_lim[i] = source_mesh.r_lim[i];
        this->theta_lim[i] = source_mesh.theta_lim[i];
        this->phi_lim[i] = source_mesh.phi_lim[i];
    }
    this->num_face = source_mesh.num_face;
    this->num_cell = source_mesh.num_cell;
    this->r_points = source_mesh.r_points;
    this->num_leaf_cells = source_mesh.num_leaf_cells;
    this->num_leaf_faces = source_mesh.num_leaf_faces;
    this->n_parameters = source_mesh.n_parameters;
    this->new_to_old_index = source_mesh.new_to_old_index;
    this->old_to_new_index = source_mesh.old_to_new_index;
    // cout<<"A"<<endl;
    // deep copy pointers
    unordered_map<Cell*, Cell*> cc_mp;  // key: old mesh, value: new mesh
    unordered_map<Face*, Face*> ff_mp;

    this->cells.resize(source_mesh.cells.size());
    for (int i = 0; i < source_mesh.cells.size(); i++) {
        this->cells[i].resize(source_mesh.cells[i].size());
        for (int j = 0; j < source_mesh.cells[i].size(); j++) {
            Cell* c = source_mesh.cells[i][j];
            this->cells[i][j] = new Cell(
                c->_r[0], c->_theta[0], c->_phi[0], c->_r[1], c->_theta[1],
                c->_phi[1], c->level, c->parameters.size(), c->isleaf);
            this->cells[i][j]->set_id(c->id);
            this->cells[i][j]->_density = c->_density;
            this->cells[i][j]->set_ordering_forward(c->get_ordering_forward());
            cc_mp[c] = this->cells[i][j];
        }
    }

    // cout<<"A"<<endl;

    this->faces.resize(source_mesh.faces.size());
    for (int i = 0; i < source_mesh.faces.size(); i++) {
        this->faces[i].resize(source_mesh.faces[i].size());
        for (int j = 0; j < source_mesh.faces[i].size(); j++) {
            Face* f = source_mesh.faces[i][j];
            this->faces[i][j] = new Face(NULL, NULL, f->rc, f->thetac, f->phic,
                                         f->direction, f->level, f->isleaf);
            ff_mp[f] = this->faces[i][j];
        }
    }

    // set associated pointers
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < cells[i].size(); j++) {
            Cell* c = source_mesh.cells[i][j];
            for (int k = 0; k < 8; k++) {
                if (c->child_cells[k] != NULL) {
                    cc_mp[c]->child_cells[k] = cc_mp[c->child_cells[k]];
                } else {
                    cc_mp[c]->child_cells[k] = NULL;
                }
            }
            // cout<<"X"<<endl;
            for (int k = 0; k < 2; k++) {
                if (c->external_faces_r[k] != NULL) {
                    cc_mp[c]->external_faces_r[k] =
                        ff_mp[c->external_faces_r[k]];
                } else {
                    cc_mp[c]->external_faces_r[k] = NULL;
                }

                if (c->external_faces_theta[k] != NULL) {
                    cc_mp[c]->external_faces_theta[k] =
                        ff_mp[c->external_faces_theta[k]];
                } else {
                    cc_mp[c]->external_faces_theta[k] = NULL;
                }

                if (c->external_faces_phi[k] != NULL) {
                    cc_mp[c]->external_faces_phi[k] =
                        ff_mp[c->external_faces_phi[k]];
                } else {
                    cc_mp[c]->external_faces_phi[k] = NULL;
                }
            }
            for (int k = 0; k < 4; k++) {
                if (c->internal_faces_r[k] != NULL) {
                    cc_mp[c]->internal_faces_r[k] =
                        ff_mp[c->internal_faces_r[k]];
                } else {
                    cc_mp[c]->internal_faces_r[k] = NULL;
                }

                if (c->internal_faces_theta[k] != NULL) {
                    cc_mp[c]->internal_faces_theta[k] =
                        ff_mp[c->internal_faces_theta[k]];
                } else {
                    cc_mp[c]->internal_faces_theta[k] = NULL;
                }

                if (c->internal_faces_phi[k] != NULL) {
                    cc_mp[c]->internal_faces_phi[k] =
                        ff_mp[c->internal_faces_phi[k]];
                } else {
                    cc_mp[c]->internal_faces_phi[k] = NULL;
                }
            }
        }
    }

    for (int i = 0; i < faces.size(); i++) {
        for (int j = 0; j < faces[i].size(); j++) {
            Face* f = source_mesh.faces[i][j];
            for (int k = 0; k < 4; k++) {
                if (f->child_faces[k] != NULL) {
                    ff_mp[f]->child_faces[k] = ff_mp[f->child_faces[k]];
                } else {
                    ff_mp[f]->child_faces[k] = NULL;
                }
            }
            for (int k = 0; k < 2; k++) {
                if (f->neigh_cells[k] != NULL) {
                    ff_mp[f]->neigh_cells[k] = cc_mp[f->neigh_cells[k]];
                } else {
                    ff_mp[f]->neigh_cells[k] = NULL;
                }
            }
        }
    }

    leaf_cells.resize(source_mesh.leaf_cells.size());
    ro_leaf_cells.resize(source_mesh.ro_leaf_cells.size());
    leaf_faces.resize(source_mesh.leaf_faces.size());
    for (int i = 0; i < leaf_cells.size(); i++) {
        Cell* c = source_mesh.leaf_cells[i];
        this->leaf_cells[i] = cc_mp[c];
    }
    for (int i = 0; i < ro_leaf_cells.size(); i++) {
        Cell* c = source_mesh.ro_leaf_cells[i];
        this->ro_leaf_cells[i] = cc_mp[c];
    }
    for (int i = 0; i < leaf_faces.size(); i++) {
        Face* f = source_mesh.leaf_faces[i];
        this->leaf_faces[i] = ff_mp[f];
    }
    return *this;
}

void Mesh::clear_all() {
    for (int i = 0; i < num_cell.size(); i++) {
        for (int j = 0; j < num_cell[i]; j++) {
            delete cells[i][j];
            cells[i][j] = NULL;
        }
    }

    for (int i = 0; i < num_face.size(); i++) {
        for (int j = 0; j < num_face[i]; j++) {
            delete faces[i][j];
            faces[i][j] = NULL;
        }
    }

    // Pointer leaf_cells[i] points to the same Cell object as
    // a certain element of vector<vector<Cell*>>cells do, so it has been
    // deleted through the above loops, and should not be deleted again.
    // So it is with leaf_faces
    for (int i = 0; i < leaf_cells.size(); i++) {
        leaf_cells[i] = NULL;
    }
    for (int j = 0; j < leaf_faces.size(); j++) {
        leaf_faces[j] = NULL;
    }

    cells.clear();
    faces.clear();
    leaf_cells.clear();
    ro_leaf_cells.clear();
    leaf_faces.clear();
    num_leaf_cells = 0;
    num_leaf_faces = 0;
    nr = 0;
    ntheta = 0;
    nphi = 0;
    num_cell.clear();
    num_face.clear();
    new_to_old_index.clear();
    old_to_new_index.clear();
}
Tesseroid& Mesh::get_elem(const unsigned int i) {
    assert(i < this->n_elems());
    assert(leaf_cells[i] != NULL);
    return *(this->leaf_cells[i]);
}
Tesseroid& Mesh::get_elem_reordered(const unsigned int i) {
    assert(i < this->n_elems());
    assert(ro_leaf_cells[i] != NULL);
    return *(this->ro_leaf_cells[i]);
}

int Mesh::get_reordered_id(const unsigned int& i) const {
    return this->new_to_old_index[i];
}
Cell* Mesh::get_element_level_0(unsigned int i_lat, unsigned int j_lon,
                                unsigned int k_r) {
    assert(i_lat < ntheta);
    assert(j_lon < nphi);
    assert(k_r < nr);
    return this->cells[0][i_lat * nphi * nr + j_lon * nr + k_r];
}

void Mesh::generate_regular_geographic_mesh(double lats[2], int n_lat,
                                            double lons[2], int n_lon,
                                            double depth_model[2], int n_depth,
                                            double reference_surface,
                                            int num_para) {
    for (int k = 0; k < 2; k++) {
        assert(great_equal(lats[k], -90) && great_equal(90.0, lats[k]));
        assert(great_equal(lons[k], -180.0) && great_equal(180.0, lons[k]));
    }
    double theta_model[2] = {90.0 - lats[1], 90.0 - lats[0]};
    double phi_model[2] = {180.0 + lons[0], 180.0 + lons[1]};
    double r_model[2] = {reference_surface - depth_model[1],
                         reference_surface - depth_model[0]};

    this->generate_regular_mesh(theta_model, n_lat, phi_model, n_lon, r_model,
                                n_depth, num_para);
}

void Mesh::generate_regular_mesh(double theta_model[2], int n_theta,
                                 double phi_model[2], int n_phi,
                                 double r_model[2], int n_r, int num_para) {
    this->clear_all();
    for (int k = 0; k < 2; k++) {
        assert(great_equal(theta_model[k], 0.0) &&
               great_equal(180.0, theta_model[k]));
        assert(great_equal(phi_model[k], 0.0) &&
               great_equal(360.0, phi_model[k]));
    }

    double theta_model_radian[2];
    double phi_model_radian[2];

    for (int i = 0; i < 2; i++) {
        theta_model_radian[i] = theta_model[i] * GS::PI / 180.0;
        phi_model_radian[i] = phi_model[i] * GS::PI / 180.0;
    }

    for (int i = 0; i < 2; i++) {
        this->r_lim[i] = r_model[i];
        this->theta_lim[i] = theta_model_radian[i];
        this->phi_lim[i] = phi_model_radian[i];
    }

    this->n_parameters = num_para;
    this->nr = n_r;
    this->ntheta = n_theta;
    this->nphi = n_phi;

    r_points = VectorXd::LinSpaced(n_r + 1, r_model[0], r_model[1]);
    VectorXd theta_points = VectorXd::LinSpaced(
        n_theta + 1, theta_model_radian[0], theta_model_radian[1]);
    VectorXd phi_points = VectorXd::LinSpaced(n_phi + 1, phi_model_radian[0],
                                              phi_model_radian[1]);

    cells.push_back(vector<Cell*>(0));
    cells[0].resize(n_r * n_theta * n_phi);

    faces.push_back(vector<Face*>(0));
    faces[0].resize(3 * n_r * n_theta * n_phi + nr * n_theta + n_theta * n_phi +
                    n_r * n_phi);

    num_cell.push_back(cells[0].size());
    num_face.push_back(faces[0].size());

    // construct cells
    int id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                cells[0][id] =
                    new Cell(r_points(k), theta_points(i), phi_points(j),
                             r_points(k + 1), theta_points(i + 1),
                             phi_points(j + 1), 0, num_para, true);
                for (int ip = 0; ip < num_para; ip++) {
                    cells[0][id]->set_parameter(0, ip);
                }
                cells[0][id]->set_id(id);
                id++;
            }
        }
    }

    num_leaf_cells = cells[0].size();
    leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        leaf_cells[i] = cells[0][i];
    }

    new_to_old_index.resize(num_leaf_cells);
    id = 0;
    assert(num_leaf_cells = (n_theta * n_phi * n_r));
    int forward = 1;  // forward: 1, inverse: 0
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                assert(id == (i * (n_phi * n_r) + j * n_r + k));
                new_to_old_index[i * (n_phi * n_r) + j * n_r + k] =
                    i * (n_phi * n_r) + j * n_r + forward * k +
                    (1 - forward) * (n_r - 1 - k);
                if (forward == 1) {
                    leaf_cells[id]->set_ordering_forward(true);
                } else {
                    leaf_cells[id]->set_ordering_forward(false);
                }
                id++;
            }
            forward = 1 - forward;
        }
    }
    old_to_new_index.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        int old_index = new_to_old_index[i];
        int new_index = i;
        old_to_new_index[old_index] = new_index;
    }
    ro_leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; ++i) {
        ro_leaf_cells[i] = leaf_cells[new_to_old_index[i]];
    }
    // construct faces

    id = 0;
    double thetac, phic, rc;
    double dtheta, dphi, dr;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r + 1; k++) {
                if (k == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc - 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (k == n_r) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + j * n_r + k - 1],
                                 NULL, RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(
                        dr, dtheta, dphi);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc + 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + j * n_r + k - 1],
                        cells[0][i * n_phi * n_r + j * n_r + k], RADIUS);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }

                if (k > 0) {
                    int index = i * n_phi * n_r + j * n_r + k - 1;
                    cells[0][index]->set_external_faces(faces[0][id - 1],
                                                        faces[0][id], RADIUS);
                }

                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta + 1; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                if (i == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], NORTH_SOUTH);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc, thetac - 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (i == n_theta) {
                    faces[0][id] =
                        new Face(cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                                 NULL, NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc, thetac + 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (i > 0) {
                    int index = (i - 1) * n_phi * n_r + j * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - n_phi * n_r], faces[0][id], NORTH_SOUTH);
                }
                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi + 1; j++) {
            for (int k = 0; k < n_r; k++) {
                if (j == 0) {
                    faces[0][id] =
                        new Face(NULL, cells[0][i * n_phi * n_r + j * n_r + k],
                                 WEST_EAST);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic - 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (j == n_phi) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                                 NULL, WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic + 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (j > 0) {
                    int index = i * n_phi * n_r + (j - 1) * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - 1 * nr], faces[0][id], WEST_EAST);
                }
                id++;
            }
        }
    }

    // for (int i = 0; i < n_x; i++)
    // {
    //     for (int j = 0; j < n_y; j++)
    //     {
    //         for (int k = 0; k < n_z; k++)
    //         {

    //         }
    //     }
    // }
    num_leaf_faces = faces[0].size();
    leaf_faces.resize(num_leaf_faces);
    for (int i = 0; i < num_leaf_faces; i++) {
        leaf_faces[i] = faces[0][i];
    }

    this->set_n_parameter(num_para);
    // this->sort(0);
}

void Mesh::generate_regular_mesh(double theta_model[2], int n_theta,
                                 double phi_model[2], int n_phi,
                                 VectorXd& r_points_, int num_para) {
    this->clear_all();
    for (int k = 0; k < 2; k++) {
        assert(great_equal(theta_model[k], 0.0) &&
               great_equal(180.0, theta_model[k]));
        assert(great_equal(phi_model[k], 0.0) &&
               great_equal(360.0, phi_model[k]));
    }

    double theta_model_radian[2];
    double phi_model_radian[2];

    for (int i = 0; i < 2; i++) {
        theta_model_radian[i] = theta_model[i] * GS::PI / 180.0;
        phi_model_radian[i] = phi_model[i] * GS::PI / 180.0;
    }

    this->r_points = r_points_;
    r_lim[0] = r_points(0);
    r_lim[1] = r_points(r_points.size() - 1);
    for (int i = 0; i < 2; i++) {
        this->theta_lim[i] = theta_model_radian[i];
        this->phi_lim[i] = phi_model_radian[i];
    }

    this->n_parameters = num_para;
    this->ntheta = n_theta;
    this->nphi = n_phi;

    int n_r = r_points.size() - 1;
    this->nr = n_r;

    VectorXd theta_points = VectorXd::LinSpaced(
        n_theta + 1, theta_model_radian[0], theta_model_radian[1]);
    VectorXd phi_points = VectorXd::LinSpaced(n_phi + 1, phi_model_radian[0],
                                              phi_model_radian[1]);

    cells.push_back(vector<Cell*>(0));
    cells[0].resize(n_r * n_theta * n_phi);

    faces.push_back(vector<Face*>(0));
    faces[0].resize(3 * n_r * n_theta * n_phi + nr * n_theta + n_theta * n_phi +
                    n_r * n_phi);

    num_cell.push_back(cells[0].size());
    num_face.push_back(faces[0].size());

    // construct cells
    int id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                cells[0][id] =
                    new Cell(r_points(k), theta_points(i), phi_points(j),
                             r_points(k + 1), theta_points(i + 1),
                             phi_points(j + 1), 0, num_para, true);
                for (int ip = 0; ip < num_para; ip++) {
                    cells[0][id]->set_parameter(0, ip);
                }
                cells[0][id]->set_id(id);
                id++;
            }
        }
    }
    num_leaf_cells = cells[0].size();
    leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        leaf_cells[i] = cells[0][i];
    }

    new_to_old_index.resize(num_leaf_cells);
    id = 0;
    assert(num_leaf_cells = (n_theta * n_phi * n_r));
    int forward = 1;  // forward: 1, inverse: 0
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                assert(id == (i * (n_phi * n_r) + j * n_r + k));
                new_to_old_index[i * (n_phi * n_r) + j * n_r + k] =
                    i * (n_phi * n_r) + j * n_r + forward * k +
                    (1 - forward) * (n_r - 1 - k);
                if (forward == 1) {
                    leaf_cells[id]->set_ordering_forward(true);
                } else {
                    leaf_cells[id]->set_ordering_forward(false);
                }
                id++;
            }
            forward = 1 - forward;
        }
    }
    old_to_new_index.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        int old_index = new_to_old_index[i];
        int new_index = i;
        old_to_new_index[old_index] = new_index;
    }
    ro_leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; ++i) {
        ro_leaf_cells[i] = leaf_cells[new_to_old_index[i]];
    }
    // construct faces

    id = 0;
    double thetac, phic, rc;
    double dtheta, dphi, dr;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r + 1; k++) {
                if (k == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc - 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (k == n_r) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + j * n_r + k - 1],
                                 NULL, RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(
                        dr, dtheta, dphi);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc + 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + j * n_r + k - 1],
                        cells[0][i * n_phi * n_r + j * n_r + k], RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc + 0.5 * dr,
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }

                if (k > 0) {
                    int index = i * n_phi * n_r + j * n_r + k - 1;
                    cells[0][index]->set_external_faces(faces[0][id - 1],
                                                        faces[0][id], RADIUS);
                }

                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta + 1; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                if (i == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], NORTH_SOUTH);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc, thetac - 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (i == n_theta) {
                    faces[0][id] =
                        new Face(cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                                 NULL, NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc, thetac + 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (i > 0) {
                    int index = (i - 1) * n_phi * n_r + j * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - n_phi * n_r], faces[0][id], NORTH_SOUTH);
                }
                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi + 1; j++) {
            for (int k = 0; k < n_r; k++) {
                if (j == 0) {
                    faces[0][id] =
                        new Face(NULL, cells[0][i * n_phi * n_r + j * n_r + k],
                                 WEST_EAST);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic - 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (j == n_phi) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                                 NULL, WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic + 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (j > 0) {
                    int index = i * n_phi * n_r + (j - 1) * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - 1 * nr], faces[0][id], WEST_EAST);
                }
                id++;
            }
        }
    }

    // for (int i = 0; i < n_x; i++)
    // {
    //     for (int j = 0; j < n_y; j++)
    //     {
    //         for (int k = 0; k < n_z; k++)
    //         {

    //         }
    //     }
    // }
    num_leaf_faces = faces[0].size();
    leaf_faces.resize(num_leaf_faces);
    for (int i = 0; i < num_leaf_faces; i++) {
        leaf_faces[i] = faces[0][i];
    }

    this->set_n_parameter(num_para);
    // this->sort(0);
}

void Mesh::generate_radially_varying_mesh_pow(
    double theta_model[2], int n_theta, double phi_model[2], int n_phi,
    double depth_1, double depth_2, int n_r, double nth_power,
    double reference_surface, int num_para) {
    this->clear_all();
    for (int k = 0; k < 2; k++) {
        assert(great_equal(theta_model[k], 0.0) &&
               great_equal(180.0, theta_model[k]));
        assert(great_equal(phi_model[k], 0.0) &&
               great_equal(360.0, phi_model[k]));
    }

    double theta_model_radian[2];
    double phi_model_radian[2];

    r_lim[0] = reference_surface - depth_2;
    r_lim[1] = reference_surface;
    for (int i = 0; i < 2; i++) {
        theta_model_radian[i] = theta_model[i] * GS::PI / 180.0;
        phi_model_radian[i] = phi_model[i] * GS::PI / 180.0;
    }

    for (int i = 0; i < 2; i++) {
        // this->r_lim[i] = r_model[i];
        this->theta_lim[i] = theta_model_radian[i];
        this->phi_lim[i] = phi_model_radian[i];
    }

    this->n_parameters = num_para;
    this->nr = n_r;
    this->ntheta = n_theta;
    this->nphi = n_phi;

    assert(!(n_r < 2));
    double x1 = pow(depth_1, 1.0 / nth_power);
    double x2 = pow(depth_2, 1.0 / nth_power);
    VectorXd x_points = VectorXd::LinSpaced(n_r, x1, x2);
    // VectorXd r_points(n_r + 1);
    this->r_points.resize(n_r + 1);
    r_points(0) = reference_surface - depth_2;
    r_points(n_r - 1) = reference_surface - depth_1;
    r_points(n_r) = reference_surface;

    for (int i = 1; i < n_r - 1; i++) {
        r_points(i) =
            reference_surface - std::pow(x_points(n_r - 1 - i), nth_power);
    }
    VectorXd theta_points = VectorXd::LinSpaced(
        n_theta + 1, theta_model_radian[0], theta_model_radian[1]);
    VectorXd phi_points = VectorXd::LinSpaced(n_phi + 1, phi_model_radian[0],
                                              phi_model_radian[1]);

    cells.push_back(vector<Cell*>(0));
    cells[0].resize(n_r * n_theta * n_phi);

    faces.push_back(vector<Face*>(0));
    faces[0].resize(3 * n_r * n_theta * n_phi + nr * n_theta + n_theta * n_phi +
                    n_r * n_phi);

    num_cell.push_back(cells[0].size());
    num_face.push_back(faces[0].size());

    // construct cells
    int id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                cells[0][id] =
                    new Cell(r_points(k), theta_points(i), phi_points(j),
                             r_points(k + 1), theta_points(i + 1),
                             phi_points(j + 1), 0, num_para, true);
                for (int ip = 0; ip < num_para; ip++) {
                    cells[0][id]->set_parameter(0, ip);
                }
                cells[0][id]->set_id(id);
                id++;
            }
        }
    }

    num_leaf_cells = cells[0].size();
    leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        leaf_cells[i] = cells[0][i];
    }

    new_to_old_index.resize(num_leaf_cells);
    id = 0;
    assert(num_leaf_cells = (n_theta * n_phi * n_r));
    int forward = 1;  // forward: 1, inverse: 0
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                assert(id == (i * (n_phi * n_r) + j * n_r + k));
                new_to_old_index[i * (n_phi * n_r) + j * n_r + k] =
                    i * (n_phi * n_r) + j * n_r + forward * k +
                    (1 - forward) * (n_r - 1 - k);
                if (forward == 1) {
                    leaf_cells[id]->set_ordering_forward(true);
                } else {
                    leaf_cells[id]->set_ordering_forward(false);
                }
                id++;
            }
            forward = 1 - forward;
        }
    }
    old_to_new_index.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        int old_index = new_to_old_index[i];
        int new_index = i;
        old_to_new_index[old_index] = new_index;
    }
    ro_leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; ++i) {
        ro_leaf_cells[i] = leaf_cells[new_to_old_index[i]];
    }
    // construct faces

    id = 0;
    double thetac, phic, rc;
    double dtheta, dphi, dr;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r + 1; k++) {
                if (k == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc - 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (k == n_r) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + j * n_r + k - 1],
                                 NULL, RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(
                        dr, dtheta, dphi);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc + 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + j * n_r + k - 1],
                        cells[0][i * n_phi * n_r + j * n_r + k], RADIUS);
                    // cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(dr,
                    // dtheta, dphi); faces[0][id]->set_center(rc + 0.5 * dr,
                    // 0.5 * (thetac + thetac1), 0.5 * (phic + phic1));
                    faces[0][id]->set_center(
                        cells[0][i * n_phi * n_r + j * n_r + k - 1]->_r[1],
                        0.5 * (thetac + thetac1), 0.5 * (phic + phic1));
                }

                if (k > 0) {
                    int index = i * n_phi * n_r + j * n_r + k - 1;
                    cells[0][index]->set_external_faces(faces[0][id - 1],
                                                        faces[0][id], RADIUS);
                }

                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta + 1; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                if (i == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], NORTH_SOUTH);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc, thetac - 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (i == n_theta) {
                    faces[0][id] =
                        new Face(cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                                 NULL, NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc, thetac + 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (i > 0) {
                    int index = (i - 1) * n_phi * n_r + j * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - n_phi * n_r], faces[0][id], NORTH_SOUTH);
                }
                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi + 1; j++) {
            for (int k = 0; k < n_r; k++) {
                if (j == 0) {
                    faces[0][id] =
                        new Face(NULL, cells[0][i * n_phi * n_r + j * n_r + k],
                                 WEST_EAST);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic - 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (j == n_phi) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                                 NULL, WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic + 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (j > 0) {
                    int index = i * n_phi * n_r + (j - 1) * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - 1 * nr], faces[0][id], WEST_EAST);
                }
                id++;
            }
        }
    }

    // for (int i = 0; i < n_x; i++)
    // {
    //     for (int j = 0; j < n_y; j++)
    //     {
    //         for (int k = 0; k < n_z; k++)
    //         {

    //         }
    //     }
    // }
    num_leaf_faces = faces[0].size();
    leaf_faces.resize(num_leaf_faces);
    for (int i = 0; i < num_leaf_faces; i++) {
        leaf_faces[i] = faces[0][i];
    }

    this->set_n_parameter(num_para);
    // this->sort(0);
}

void Mesh::generate_radially_varying_mesh(double theta_model[2], int n_theta,
                                          double phi_model[2], int n_phi,
                                          double h1, double depth, int n_r,
                                          double reference_surface,
                                          int num_para) {
    this->clear_all();
    for (int k = 0; k < 2; k++) {
        assert(great_equal(theta_model[k], 0.0) &&
               great_equal(180.0, theta_model[k]));
        assert(great_equal(phi_model[k], 0.0) &&
               great_equal(360.0, phi_model[k]));
    }

    double theta_model_radian[2];
    double phi_model_radian[2];

    r_lim[0] = reference_surface - depth;
    r_lim[1] = reference_surface;
    for (int i = 0; i < 2; i++) {
        theta_model_radian[i] = theta_model[i] * GS::PI / 180.0;
        phi_model_radian[i] = phi_model[i] * GS::PI / 180.0;
    }

    for (int i = 0; i < 2; i++) {
        // this->r_lim[i] = r_model[i];
        this->theta_lim[i] = theta_model_radian[i];
        this->phi_lim[i] = phi_model_radian[i];
    }

    this->n_parameters = num_para;
    this->nr = n_r;
    this->ntheta = n_theta;
    this->nphi = n_phi;

    assert(!(n_r < 2));

    double hn = (2.0 * depth) / n_r - h1;
    double dh = (hn - h1) / (n_r - 1);
    // VectorXd r_points(n_r + 1);
    this->r_points.resize(n_r + 1);
    r_points(0) = reference_surface - depth;
    r_points(n_r) = reference_surface;

    cout << (h1 + (n_r - 1) * dh) << endl;
    for (int i = 1; i < n_r; i++) {
        double depth_i = (h1 + (h1 + (n_r - i - 1) * dh)) * (n_r - i) / 2;
        cout << (h1 + (n_r - i - 1) * dh) << endl;
        r_points(i) = reference_surface - depth_i;
    }
    cout << r_points(0) << ", " << r_points(n_r) << endl;
    VectorXd theta_points = VectorXd::LinSpaced(
        n_theta + 1, theta_model_radian[0], theta_model_radian[1]);
    VectorXd phi_points = VectorXd::LinSpaced(n_phi + 1, phi_model_radian[0],
                                              phi_model_radian[1]);

    cells.push_back(vector<Cell*>(0));
    cells[0].resize(n_r * n_theta * n_phi);

    faces.push_back(vector<Face*>(0));
    faces[0].resize(3 * n_r * n_theta * n_phi + nr * n_theta + n_theta * n_phi +
                    n_r * n_phi);

    num_cell.push_back(cells[0].size());
    num_face.push_back(faces[0].size());

    // construct cells
    int id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                cells[0][id] =
                    new Cell(r_points(k), theta_points(i), phi_points(j),
                             r_points(k + 1), theta_points(i + 1),
                             phi_points(j + 1), 0, num_para, true);
                for (int ip = 0; ip < num_para; ip++) {
                    cells[0][id]->set_parameter(0, ip);
                }
                cells[0][id]->set_id(id);
                id++;
            }
        }
    }

    num_leaf_cells = cells[0].size();
    leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        leaf_cells[i] = cells[0][i];
    }

    new_to_old_index.resize(num_leaf_cells);
    id = 0;
    assert(num_leaf_cells = (n_theta * n_phi * n_r));
    int forward = 1;  // forward: 1, inverse: 0
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                assert(id == (i * (n_phi * n_r) + j * n_r + k));
                new_to_old_index[i * (n_phi * n_r) + j * n_r + k] =
                    i * (n_phi * n_r) + j * n_r + forward * k +
                    (1 - forward) * (n_r - 1 - k);
                if (forward == 1) {
                    leaf_cells[id]->set_ordering_forward(true);
                } else {
                    leaf_cells[id]->set_ordering_forward(false);
                }
                id++;
            }
            forward = 1 - forward;
        }
    }
    old_to_new_index.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        int old_index = new_to_old_index[i];
        int new_index = i;
        old_to_new_index[old_index] = new_index;
    }
    ro_leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; ++i) {
        ro_leaf_cells[i] = leaf_cells[new_to_old_index[i]];
    }
    // construct faces

    id = 0;
    double thetac, phic, rc;
    double dtheta, dphi, dr;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r + 1; k++) {
                if (k == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc - 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (k == n_r) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + j * n_r + k - 1],
                                 NULL, RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(
                        dr, dtheta, dphi);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc + 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    // cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(dr,
                    // dtheta, dphi);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + j * n_r + k - 1],
                        cells[0][i * n_phi * n_r + j * n_r + k], RADIUS);
                    // double f_rc=
                    faces[0][id]->set_center(
                        cells[0][i * n_phi * n_r + j * n_r + k - 1]->_r[1],
                        0.5 * (thetac + thetac1), 0.5 * (phic + phic1));
                }

                if (k > 0) {
                    int index = i * n_phi * n_r + j * n_r + k - 1;
                    cells[0][index]->set_external_faces(faces[0][id - 1],
                                                        faces[0][id], RADIUS);
                }

                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta + 1; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                if (i == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], NORTH_SOUTH);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc, thetac - 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (i == n_theta) {
                    faces[0][id] =
                        new Face(cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                                 NULL, NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc, thetac + 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (i > 0) {
                    int index = (i - 1) * n_phi * n_r + j * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - n_phi * n_r], faces[0][id], NORTH_SOUTH);
                }
                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi + 1; j++) {
            for (int k = 0; k < n_r; k++) {
                if (j == 0) {
                    faces[0][id] =
                        new Face(NULL, cells[0][i * n_phi * n_r + j * n_r + k],
                                 WEST_EAST);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic - 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (j == n_phi) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                                 NULL, WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic + 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (j > 0) {
                    int index = i * n_phi * n_r + (j - 1) * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - 1 * nr], faces[0][id], WEST_EAST);
                }
                id++;
            }
        }
    }

    // for (int i = 0; i < n_x; i++)
    // {
    //     for (int j = 0; j < n_y; j++)
    //     {
    //         for (int k = 0; k < n_z; k++)
    //         {

    //         }
    //     }
    // }
    num_leaf_faces = faces[0].size();
    leaf_faces.resize(num_leaf_faces);
    for (int i = 0; i < num_leaf_faces; i++) {
        leaf_faces[i] = faces[0][i];
    }

    this->set_n_parameter(num_para);
    // this->sort(0);
}

void Mesh::generate_radially_varying_mesh_log(
    double theta_model[2], int n_theta, double phi_model[2], int n_phi,
    double depth_1, double depth_2, int n_r, double base,
    double reference_surface, int num_para) {
    this->clear_all();
    for (int k = 0; k < 2; k++) {
        assert(great_equal(theta_model[k], 0.0) &&
               great_equal(180.0, theta_model[k]));
        assert(great_equal(phi_model[k], 0.0) &&
               great_equal(360.0, phi_model[k]));
    }

    double theta_model_radian[2];
    double phi_model_radian[2];

    r_lim[0] = reference_surface - depth_2;
    r_lim[1] = reference_surface;
    for (int i = 0; i < 2; i++) {
        theta_model_radian[i] = theta_model[i] * GS::PI / 180.0;
        phi_model_radian[i] = phi_model[i] * GS::PI / 180.0;
    }

    for (int i = 0; i < 2; i++) {
        // this->r_lim[i] = r_model[i];
        this->theta_lim[i] = theta_model_radian[i];
        this->phi_lim[i] = phi_model_radian[i];
    }

    this->n_parameters = num_para;
    this->nr = n_r;
    this->ntheta = n_theta;
    this->nphi = n_phi;

    assert(!(n_r < 2));
    double a = log(depth_1) / log(base);
    double b = log(depth_2) / log(base);
    double log_interval = (b - a) / (n_r - 1);

    // VectorXd r_points(n_r + 1);
    this->r_points.resize(n_r + 1);

    r_points(0) = reference_surface - depth_2;
    r_points(n_r) = reference_surface;
    for (int i = 1; i < nr; i++) {
        double depth = pow(base, a + (n_r - i - 1) * log_interval);
        r_points(i) = reference_surface - depth;
    }
    VectorXd theta_points = VectorXd::LinSpaced(
        n_theta + 1, theta_model_radian[0], theta_model_radian[1]);
    VectorXd phi_points = VectorXd::LinSpaced(n_phi + 1, phi_model_radian[0],
                                              phi_model_radian[1]);

    cells.push_back(vector<Cell*>(0));
    cells[0].resize(n_r * n_theta * n_phi);

    faces.push_back(vector<Face*>(0));
    faces[0].resize(3 * n_r * n_theta * n_phi + nr * n_theta + n_theta * n_phi +
                    n_r * n_phi);

    num_cell.push_back(cells[0].size());
    num_face.push_back(faces[0].size());

    // construct cells
    int id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                cells[0][id] =
                    new Cell(r_points(k), theta_points(i), phi_points(j),
                             r_points(k + 1), theta_points(i + 1),
                             phi_points(j + 1), 0, num_para, true);
                for (int ip = 0; ip < num_para; ip++) {
                    cells[0][id]->set_parameter(0, ip);
                }
                cells[0][id]->set_id(id);
                id++;
            }
        }
    }
    num_leaf_cells = cells[0].size();
    leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        leaf_cells[i] = cells[0][i];
    }
    new_to_old_index.resize(num_leaf_cells);
    id = 0;
    assert(num_leaf_cells = (n_theta * n_phi * n_r));
    int forward = 1;  // forward: 1, inverse: 0
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                assert(id == (i * (n_phi * n_r) + j * n_r + k));
                new_to_old_index[i * (n_phi * n_r) + j * n_r + k] =
                    i * (n_phi * n_r) + j * n_r + forward * k +
                    (1 - forward) * (n_r - 1 - k);
                if (forward == 1) {
                    leaf_cells[id]->set_ordering_forward(true);
                } else {
                    leaf_cells[id]->set_ordering_forward(false);
                }
                id++;
            }
            forward = 1 - forward;
        }
    }
    old_to_new_index.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; i++) {
        int old_index = new_to_old_index[i];
        int new_index = i;
        old_to_new_index[old_index] = new_index;
    }
    ro_leaf_cells.resize(num_leaf_cells);
    for (int i = 0; i < num_leaf_cells; ++i) {
        ro_leaf_cells[i] = leaf_cells[new_to_old_index[i]];
    }
    // construct faces

    id = 0;
    double thetac, phic, rc;
    double dtheta, dphi, dr;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r + 1; k++) {
                if (k == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc - 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (k == n_r) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + j * n_r + k - 1],
                                 NULL, RADIUS);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(
                        dr, dtheta, dphi);
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc + 0.5 * dr, thetac, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + j * n_r + k - 1],
                        cells[0][i * n_phi * n_r + j * n_r + k], RADIUS);
                    // cells[0][i * n_phi * n_r + j * n_r + k - 1]->get_size(dr,
                    // dtheta, dphi); faces[0][id]->set_center(rc + 0.5 * dr,
                    // 0.5 * (thetac + thetac1), 0.5 * (phic + phic1));
                    faces[0][id]->set_center(
                        cells[0][i * n_phi * n_r + j * n_r + k - 1]->_r[1],
                        0.5 * (thetac + thetac1), 0.5 * (phic + phic1));
                }

                if (k > 0) {
                    int index = i * n_phi * n_r + j * n_r + k - 1;
                    cells[0][index]->set_external_faces(faces[0][id - 1],
                                                        faces[0][id], RADIUS);
                }

                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta + 1; i++) {
        for (int j = 0; j < n_phi; j++) {
            for (int k = 0; k < n_r; k++) {
                if (i == 0) {
                    faces[0][id] = new Face(
                        NULL, cells[0][i * n_phi * n_r + j * n_r], NORTH_SOUTH);
                    cells[0][i * n_phi * n_r + j * n_r]->get_size(dr, dtheta,
                                                                  dphi);
                    cells[0][i * n_phi * n_r + j * n_r]->get_center(rc, thetac,
                                                                    phic);
                    faces[0][id]->set_center(rc, thetac - 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (i == n_theta) {
                    faces[0][id] =
                        new Face(cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                                 NULL, NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    faces[0][id]->set_center(rc, thetac + 0.5 * dtheta, phic);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][(i - 1) * n_phi * n_r + j * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], NORTH_SOUTH);
                    cells[0][(i - 1) * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (i > 0) {
                    int index = (i - 1) * n_phi * n_r + j * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - n_phi * n_r], faces[0][id], NORTH_SOUTH);
                }
                id++;
            }
        }
    }

    // id = 0;
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi + 1; j++) {
            for (int k = 0; k < n_r; k++) {
                if (j == 0) {
                    faces[0][id] =
                        new Face(NULL, cells[0][i * n_phi * n_r + j * n_r + k],
                                 WEST_EAST);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic - 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] == NULL &&
                           faces[0][id]->neigh_cells[1] != NULL);
                } else if (j == n_phi) {
                    faces[0][id] =
                        new Face(cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                                 NULL, WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_size(
                        dr, dtheta, dphi);
                    faces[0][id]->set_center(rc, thetac, phic + 0.5 * dphi);
                    assert(faces[0][id]->neigh_cells[0] != NULL &&
                           faces[0][id]->neigh_cells[1] == NULL);
                } else {
                    double rc1, thetac1, phic1;
                    faces[0][id] = new Face(
                        cells[0][i * n_phi * n_r + (j - 1) * n_r + k],
                        cells[0][i * n_phi * n_r + j * n_r + k], WEST_EAST);
                    cells[0][i * n_phi * n_r + (j - 1) * n_r + k]->get_center(
                        rc, thetac, phic);
                    cells[0][i * n_phi * n_r + j * n_r + k]->get_center(
                        rc1, thetac1, phic1);
                    faces[0][id]->set_center(0.5 * (rc + rc1),
                                             0.5 * (thetac + thetac1),
                                             0.5 * (phic + phic1));
                }
                if (j > 0) {
                    int index = i * n_phi * n_r + (j - 1) * n_r + k;
                    cells[0][index]->set_external_faces(
                        faces[0][id - 1 * nr], faces[0][id], WEST_EAST);
                }
                id++;
            }
        }
    }

    // for (int i = 0; i < n_x; i++)
    // {
    //     for (int j = 0; j < n_y; j++)
    //     {
    //         for (int k = 0; k < n_z; k++)
    //         {

    //         }
    //     }
    // }
    num_leaf_faces = faces[0].size();
    leaf_faces.resize(num_leaf_faces);
    for (int i = 0; i < num_leaf_faces; i++) {
        leaf_faces[i] = faces[0][i];
    }

    this->set_n_parameter(num_para);
    // this->sort(0);
}
void Mesh::set_parameter_in_a_region(double lat0, double lat1, double lon0,
                                     double lon1, double d0, double d1,
                                     double para_value,
                                     double reference_surface, int i_th) {
    for (int i = 0; i < cells[0].size(); i++) {
        Cell* c = cells[0][i];
        double cd0 = reference_surface - c->_r[1];
        double cd1 = reference_surface - c->_r[0];
        double clat0 = 90.0 - (c->_theta[1]) * 180.0 / GS::PI;
        double clat1 = 90.0 - (c->_theta[0]) * 180.0 / GS::PI;
        double clon0 = -180.0 + (c->_phi[0]) * 180.0 / GS::PI;
        double clon1 = -180.0 + (c->_phi[1]) * 180.0 / GS::PI;

        if ((cd0 > d0 || abs(cd0 - d0) < 1e-10) &&
            (cd1 < d1 || abs(cd1 - d1) < 1e-10)) {
            if ((clat0 > lat0 || abs(clat0 - lat0) < 1e-10) &&
                (clat1 < lat1 || abs(clat1 - lat1) < 1e-10)) {
                if ((clon0 > lon0 || abs(clon0 - lon0) < 1e-10) &&
                    (clon1 < lon1 || abs(clon1 - lon1) < 1e-10)) {
                    c->set_parameter(para_value, i_th);
                }
            }
        }
    }
}

void Mesh::set_block_parameter(unsigned int i_min, unsigned int i_max,
                               unsigned int j_min, unsigned int j_max,
                               unsigned int k_min, unsigned int k_max,
                               double para_value, int i_th) {
    assert(i_min >= 0 && i_max < ntheta);
    assert(j_min >= 0 && j_max < nphi);
    assert(k_min >= 0 && k_max < nr);
    for (int i = i_min; i <= i_max; i++) {
        for (int j = j_min; j <= j_max; j++) {
            for (int k = k_min; k <= k_max; k++) {
                int id = i * nphi * nr + j * nr + k;
                cells[0][id]->set_parameter(para_value, i_th);
            }
        }
    }
}

void Mesh::convert_txt_2_vtk(string input_txtfile, string output_vtkfile,
                             int n) {
    ifstream in_txt(input_txtfile);
    ofstream out_vtk(output_vtkfile);

    assert(in_txt.good());

    string line;
    // skip two lines
    for (int i = 0; i < n; i++) {
        std::getline(in_txt, line);
    }

    Mesh txt_mesh;
    txt_mesh.cells.resize(1);
    while (std::getline(in_txt, line)) {
        std::istringstream iss(line);
        double longitude0, longitude1, latitude0, latitude1, radius0, radius1;
        double temp;
        double value;
        iss >> longitude0 >> longitude1 >> latitude0 >> latitude1 >> radius0 >>
            radius1 >> temp >> temp >> temp >> temp >> value;
        double r0, theta0, phi0;
        double r1, theta1, phi1;
        r0 = radius0 * 1000.0;
        r1 = radius1 * 1000.0;
        phi0 = (longitude0 + 180.0) * GS::PI / 180.0;
        phi1 = (longitude1 + 180.0) * GS::PI / 180.0;
        theta0 = (90.0 - latitude1) * GS::PI / 180.0;
        theta1 = (90.0 - latitude0) * GS::PI / 180.0;
        Cell* c = new Cell(r0, theta0, phi0, r1, theta1, phi1);
        c->set_parameter(value);
        txt_mesh.cells[0].push_back(c);
    }
    txt_mesh.num_cell.push_back(txt_mesh.cells[0].size());
    txt_mesh.num_leaf_cells = txt_mesh.cells[0].size();
    cout << txt_mesh.num_leaf_cells << endl;
    txt_mesh.leaf_cells.resize(txt_mesh.num_leaf_cells);
    for (int i = 0; i < txt_mesh.num_leaf_cells; i++) {
        txt_mesh.leaf_cells[i] = txt_mesh.cells[0][i];
    }
    txt_mesh.out_model_vtk(output_vtkfile);
    // string out_vtkfile_linear_projection = output_vtkfile;
    // auto pos = output_vtkfile.find_first_of('.');
    // out_vtkfile_linear_projection.insert(pos, "_linear_projection");
    // txt_mesh.out_model_vtk_linear_projection(out_vtkfile_linear_projection);
}
void Mesh::out_model_txt(string filename, int ith_para,
                         double reference_surface) {
    ofstream outfile(filename.c_str());
    outfile << "#The first 6 parameters give the dimensions of a cell. The 7-9 "
               "parameters are the cell center. The 10th parameter is the "
               "value within the cell."
            << endl;
    outfile << setw(23) << left << "#Longitude 0(degree)" << setw(23) << left
            << "Longitude 1(degree)" << setw(23) << left
            << "#Latitude 0(degree)" << setw(23) << left << "Latitude 1(degree)"
            << setw(23) << left << "#R 0(km)" << setw(23) << left << "R 1(km)"
            << setw(23) << left << "Longitude c(degree)" << setw(23)
            << "Latitude c(degree)" << setw(23) << "R c(km)";
    outfile << setw(23) << left << "Depth (km)";
    outfile << setw(23) << left << "value" << endl;
    outfile << scientific;
    for (int i = 0; i < n_elems(); i++) {
        double longi1, lat1, r1;
        double longi2, lat2, r2;
        double longic, latc, rc, depthc;
        double value;

        longi1 = (this->get_elem(i)._phi[0] * 180.0 / GS::PI) - 180.0;
        longi2 = (this->get_elem(i)._phi[1] * 180.0 / GS::PI) - 180.0;
        lat1 = 90 - 180 * this->get_elem(i)._theta[1] / GS::PI;
        lat2 = 90 - 180 * this->get_elem(i)._theta[0] / GS::PI;
        r1 = this->get_elem(i)._r[0] / 1000.0;
        r2 = this->get_elem(i)._r[1] / 1000.0;

        longic = 0.5 * (longi1 + longi2);
        latc = 0.5 * (lat1 + lat2);
        rc = 0.5 * (r1 + r2);
        double depth = 0;
        depth = reference_surface / 1000.0 - rc;

        value = leaf_cells[i]->get_parameter(ith_para);
        outfile << setw(23) << left << setprecision(15) << longi1;
        outfile << setw(23) << left << setprecision(15) << longi2;
        outfile << setw(23) << left << setprecision(15) << lat1;
        outfile << setw(23) << left << setprecision(15) << lat2;
        outfile << setw(23) << left << setprecision(15) << r1;
        outfile << setw(23) << left << setprecision(15) << r2;
        outfile << setw(23) << left << setprecision(15) << longic;
        outfile << setw(23) << left << setprecision(15) << latc;
        outfile << setw(23) << left << setprecision(15) << rc;
        outfile << setw(23) << left << setprecision(15) << depth;
        outfile << setw(23) << left << setprecision(15) << value;
        outfile << endl;
    }
    cout << "The model has been written to text file: " << filename << endl;
}

#pragma optimize("", off)
void Mesh::out_model_vtk(string filename, int n,
                         vector<string> parameter_name) {
    // prepare the node
    std::set<Point> v_set;
    std::vector<std::vector<Point>> quad_prism(this->n_elems());
    // for (int i = 0; i < n_elems() && (this->get_elem(i)._phi[0] > 0.); i++)
    // for (int i = 0; i < n_elems(); i++) {
    //    if (!great_equal(this->get_elem(i)._phi[0], 0.)) {
    //        cout << get_elem(i) << endl;
    //        abort();
    //    }
    //}
    for (int i = 0; i < n_elems() && great_equal(this->get_elem(i)._phi[0], 0.);
         i++) {
        // build the point (8 points)
        Point v[20];
        double r1, r2, theta_1, theta_2, phi_1, phi_2;
        r1 = this->get_elem(i)._r[0];
        r2 = this->get_elem(i)._r[1];
        theta_1 = this->get_elem(i)._theta[0];
        theta_2 = this->get_elem(i)._theta[1];
        phi_1 = this->get_elem(i)._phi[0];
        phi_2 = this->get_elem(i)._phi[1];

        v[0] = Point(r1, theta_1, phi_1);
        v[1] = Point(r1, theta_2, phi_1);
        v[2] = Point(r1, theta_2, phi_2);
        v[3] = Point(r1, theta_1, phi_2);
        v[4] = Point(r2, theta_1, phi_1);
        v[5] = Point(r2, theta_2, phi_1);
        v[6] = Point(r2, theta_2, phi_2);
        v[7] = Point(r2, theta_1, phi_2);
        v[8] = (v[0] + v[1]) * 0.5;
        v[9] = (v[1] + v[2]) * 0.5;
        v[10] = (v[2] + v[3]) * 0.5;
        v[11] = (v[3] + v[0]) * 0.5;
        v[12] = v[8];
        v[12](0) = r2;
        v[13] = v[9];
        v[13](0) = r2;
        v[14] = v[10];
        v[14](0) = r2;
        v[15] = v[11];
        v[15](0) = r2;
        v[16] = (v[0] + v[4]) * 0.5;
        v[17] = (v[1] + v[5]) * 0.5;
        v[18] = (v[2] + v[6]) * 0.5;
        v[19] = (v[3] + v[7]) * 0.5;
        for (int j = 0; j < 20; j++) {
            v_set.insert(v[j]);
        }
        quad_prism[i].clear();
        for (int j = 0; j < 20; j++) {
            quad_prism[i].push_back(v[j]);
        }
    }
    std::map<Point, unsigned int> v_id_map;
    unsigned int counter = 0;
    for (std::set<Point>::iterator it = v_set.begin(); it != v_set.end();
         it++) {
        v_id_map[(*it)] = counter;
        counter++;
    }
    const unsigned int total_points = v_set.size();
    const unsigned int total_cells = this->n_elems();
    //----------------------------------------------------------------------
    // Open the of stream
    std::ofstream vtk_mesh(filename.c_str());
    if (!vtk_mesh.good()) {
        std::cerr << "Can not open file:\t" << filename + ".vtk" << std::endl;
    } else {
        // Parts 1-2-3, mandatory
        vtk_mesh
            << "# vtk DataFile Version 3.0\n"  // File version and identifier
            // Header info, doublely cool data
            << "Gravity inversion using tesseroids\n"
            << "ASCII\n";  // ASCII data (not BINARY)

        // Part 4, Geometry/topology, unstructured mesh
        vtk_mesh << "DATASET UNSTRUCTURED_GRID\n";  // topography and geometry

        // POINTS info (0-->n-1)
        vtk_mesh << "\nPOINTS\t" << total_points << "\tdouble\n";
        // Loop POINTS to write out coordinates
        for (std::set<Point>::iterator it = v_set.begin(); it != v_set.end();
             it++) {
            long double x, y, z;
            x = y = z = 0;
            double r = (*it)(0);
            double theta = (*it)(1);
            double phi = (*it)(2);
            rthetaphi_xyz(r, theta, phi, x, y, z);
            vtk_mesh << x << "\t"   // x-coordinate
                     << y << "\t"   // y-coordinate
                     << z << "\n";  // z-coordinate
        }

        // CELL info (0-->m-1)
        typedef std::map<Point, unsigned int>::iterator IT;

        vtk_mesh << "\nCELLS\t" << total_cells << "\t" << total_cells * (20 + 1)
                 << "\n";
        // for (unsigned int i = 0; i < total_cells &&
        // (this->get_elem(i)._phi[0] > 0.); i++)
        for (unsigned int i = 0;
             i < total_cells && great_equal(this->get_elem(i)._phi[0], 0.);
             i++) {
            if (!great_equal(this->get_elem(i)._phi[0], 0.)) {
                cout << get_elem(i) << endl;
            }
            std::vector<Point>& T = quad_prism[i];  // 20 vertex
            assert(T.size() == 20);
            unsigned int T_ID[20];
            for (int j = 0; j < 20; j++) {
                IT it = v_id_map.find(T[j]);
                assert(it != v_id_map.end());
                T_ID[j] = (*it).second;
            }
            vtk_mesh << (unsigned int)20 << "\t";
            for (int j = 0; j < 20; j++) vtk_mesh << T_ID[j] << "\t";
            vtk_mesh << "\n";
        }

        // CELL types (m)
        vtk_mesh << "\nCELL_TYPES\t" << total_cells << "\n";
        for (unsigned int i = 0; i < total_cells; i++) {
            vtk_mesh << (unsigned int)25 << "\n";  // 25-Quadratic_hexahedron
            // figure 3 in vtk format file.
        }

        // Part 5, attributes
        vtk_mesh << "\nCELL_DATA\t" << total_cells << "\n";
        for (int j = 0; j < n; j++) {
            vtk_mesh << "SCALARS " << parameter_name[j] << " double 1\n"
                     << "LOOKUP_TABLE "
                     << "table" << j << endl;
            for (unsigned int i = 0; i < total_cells; i++) {
                double value = leaf_cells[i]->get_parameter(j);
                vtk_mesh << value << "\n";
            }
        }

        /*
           vtk_mesh<<"SCALARS Elevation double 1\n"
           <<"LOOKUP_TABLE default\n";
           for(unsigned int i=0; i<total_cells; i++) {
           double color_value = this->get_elem(i)._r[1]-GS::MER;
           vtk_mesh<<color_value <<"\n";
           }
           */
        vtk_mesh << "\n";

    }  // file opened successfully

    vtk_mesh.close();
    cout << "The model has been written to vtk file: " << filename << endl;
}
#pragma optimize("", on)

#pragma optimize("", off)
void Mesh::out_model_vtk_linear_projection(string filename, int n,
                                           double reference_surface,
                                           vector<string> parameter_name) {
    // prepare the node
    std::set<Point> v_set;
    std::vector<std::vector<Point>> prism(this->n_elems());
    // for (int i = 0; i < n_elems() && (this->get_elem(i)._phi[0] > 0.); i++)
    // for (int i = 0; i < n_elems(); i++) {
    //    if (!great_equal(this->get_elem(i)._phi[0], 0.)) {
    //        cout << get_elem(i) << endl;
    //        abort();
    //    }
    //}
    for (int i = 0; i < n_elems() && great_equal(this->get_elem(i)._phi[0], 0.);
         i++) {
        // build the point (8 points)
        Point v[8];
        double r1, r2, theta_1, theta_2, phi_1, phi_2;
        r1 = this->get_elem(i)._r[0];
        r2 = this->get_elem(i)._r[1];
        theta_1 = this->get_elem(i)._theta[0];
        theta_2 = this->get_elem(i)._theta[1];
        phi_1 = this->get_elem(i)._phi[0];
        phi_2 = this->get_elem(i)._phi[1];

        v[0] = Point(r1, theta_1, phi_1);
        v[1] = Point(r1, theta_2, phi_1);
        v[2] = Point(r1, theta_2, phi_2);
        v[3] = Point(r1, theta_1, phi_2);
        v[4] = Point(r2, theta_1, phi_1);
        v[5] = Point(r2, theta_2, phi_1);
        v[6] = Point(r2, theta_2, phi_2);
        v[7] = Point(r2, theta_1, phi_2);
        for (int j = 0; j < 8; j++) {
            v_set.insert(v[j]);
        }
        prism[i].clear();
        for (int j = 0; j < 8; j++) {
            prism[i].push_back(v[j]);
        }
    }
    std::map<Point, unsigned int> v_id_map;
    unsigned int counter = 0;
    for (std::set<Point>::iterator it = v_set.begin(); it != v_set.end();
         it++) {
        v_id_map[(*it)] = counter;
        counter++;
    }
    const unsigned int total_points = v_set.size();
    const unsigned int total_cells = this->n_elems();
    //----------------------------------------------------------------------
    // Open the of stream
    std::ofstream vtk_mesh(filename.c_str());
    if (!vtk_mesh.good()) {
        std::cerr << "Can not open file:\t" << filename + ".vtk" << std::endl;
    } else {
        // Parts 1-2-3, mandatory
        vtk_mesh
            << "# vtk DataFile Version 3.0\n"  // File version and identifier
            // Header info, doublely cool data
            << "Gravity inversion using tesseroids in spherical coordinate. "
               "For convenience of slicing the model, the result is written in "
               "this file using linear projection for latitudes and "
               "longitudes.\n"
            << "ASCII\n";  // ASCII data (not BINARY)

        // Part 4, Geometry/topology, unstructured mesh
        vtk_mesh << "DATASET UNSTRUCTURED_GRID\n";  // topography and geometry

        // POINTS info (0-->n-1)
        vtk_mesh << "\nPOINTS\t" << total_points << "\tdouble\n";
        // Loop POINTS to write out coordinates
        for (std::set<Point>::iterator it = v_set.begin(); it != v_set.end();
             it++) {
            long double x, y, z;
            x = y = z = 0;
            double r = (*it)(0);
            double theta = (*it)(1);
            double phi = (*it)(2);
            x = phi * 180.0 / GS::PI - 180;
            y = 90.0 - theta * 180 / GS::PI;
            z = (r - reference_surface) / 1000.0;  // km
            vtk_mesh << x << "\t"                  // x-coordinate
                     << y << "\t"                  // y-coordinate
                     << z << "\n";                 // z-coordinate
        }

        // CELL info (0-->m-1)
        typedef std::map<Point, unsigned int>::iterator IT;

        vtk_mesh << "\nCELLS\t" << total_cells << "\t" << total_cells * (8 + 1)
                 << "\n";
        // for (unsigned int i = 0; i < total_cells &&
        // (this->get_elem(i)._phi[0] > 0.); i++)
        for (unsigned int i = 0;
             i < total_cells && great_equal(this->get_elem(i)._phi[0], 0.);
             i++) {
            if (!great_equal(this->get_elem(i)._phi[0], 0.)) {
                cout << get_elem(i) << endl;
            }
            std::vector<Point>& T = prism[i];  // 8 vertex
            assert(T.size() == 8);
            unsigned int T_ID[8];
            for (int j = 0; j < 8; j++) {
                IT it = v_id_map.find(T[j]);
                assert(it != v_id_map.end());
                T_ID[j] = (*it).second;
            }
            vtk_mesh << (unsigned int)8 << "\t";
            for (int j = 0; j < 8; j++) vtk_mesh << T_ID[j] << "\t";
            vtk_mesh << "\n";
        }

        // CELL types (m)
        vtk_mesh << "\nCELL_TYPES\t" << total_cells << "\n";
        for (unsigned int i = 0; i < total_cells; i++) {
            vtk_mesh << (unsigned int)12 << "\n";  // 12-hexahedron
            // figure 3 in vtk format file.
        }

        // Part 5, attributes
        vtk_mesh << "\nCELL_DATA\t" << total_cells << "\n";
        for (int j = 0; j < n; j++) {
            vtk_mesh << "SCALARS " << parameter_name[j] << " double 1\n"
                     << "LOOKUP_TABLE "
                     << "table" << j << endl;
            for (unsigned int i = 0; i < total_cells; i++) {
                double value = leaf_cells[i]->get_parameter(j);
                vtk_mesh << value << "\n";
            }
        }

        /*
           vtk_mesh<<"SCALARS Elevation double 1\n"
           <<"LOOKUP_TABLE default\n";
           for(unsigned int i=0; i<total_cells; i++) {
           double color_value = this->get_elem(i)._r[1]-GS::MER;
           vtk_mesh<<color_value <<"\n";
           }
           */
        vtk_mesh << "\n";

    }  // file opened successfully

    vtk_mesh.close();
    cout << "The model has been written to vtk file: " << filename << endl;
}
#pragma optimize("", on)

map<unsigned int, Cell*> Mesh::refinement(Cell* c) {
    // cout<<"count="<<std::count(leaf_cells.begin(), leaf_cells.end(),
    // c)<<endl; cout<<"isleaf="<<c->isleaf<<endl;
    assert(std::count(leaf_cells.begin(), leaf_cells.end(), c) == 1);
    bool flag = false;
    bool flag1[2] = {true, true};
    bool flag2[2] = {true, true};
    bool flag3[2] = {true, true};
    int nei_index[2] = {0, 1};
    map<unsigned int, Cell*> split_cells;
    // cout << "AB" << endl;
    do {
        for (int i = 0; i < 2; i++) {
            Face* fx = c->external_faces_theta[i];
            if ((fx->neigh_cells[0] != NULL) && (fx->neigh_cells[1] != NULL)) {
                // flag = flag && ((fx->isleaf == false) ||
                // (fx->neigh_cells[nei_index[i]]->level >= c->level));
                flag1[i] = (fx->isleaf == false) ||
                           (fx->neigh_cells[nei_index[i]]->level >= c->level);
            }

            Face* fy = c->external_faces_phi[i];
            if ((fy->neigh_cells[0] != NULL) && (fy->neigh_cells[1] != NULL)) {
                // flag = flag && ((fy->isleaf == false) ||
                // fy->neigh_cells[nei_index[i]]->level >= c->level);
                flag2[i] = (fy->isleaf == false) ||
                           fy->neigh_cells[nei_index[i]]->level >= c->level;
            }

            Face* fz = c->external_faces_r[i];
            if ((fz->neigh_cells[0] != NULL) && (fz->neigh_cells[1] != NULL)) {
                flag3[i] = (fz->isleaf == false) ||
                           fz->neigh_cells[nei_index[i]]->level >= c->level;
            }
        }
        flag = flag1[0] && flag1[1] && flag2[0] && flag2[1] && flag3[0] &&
               flag3[1];

        if (flag == false) {
            // refine neibouring cells with lower level
            for (int i = 0; i < 2; i++) {
                // cout<<i<<endl;
                Face* fx = c->external_faces_theta[i];
                if (flag1[i] == false) {
                    if (fx->neigh_cells[i] != c && fx->neigh_cells[i] != NULL) {
                        map<unsigned int, Cell*> temp;
                        temp = refinement(fx->neigh_cells[i]);
                        split_cells.insert(temp.begin(), temp.end());
                    }
                    assert(fx->neigh_cells[i]->level >= c->level);
                }

                Face* fy = c->external_faces_phi[i];
                if (flag2[i] == false) {
                    if (fy->neigh_cells[i] != c && fy->neigh_cells[i] != NULL) {
                        map<unsigned int, Cell*> temp;
                        temp = refinement(fy->neigh_cells[i]);
                        split_cells.insert(temp.begin(), temp.end());
                    }
                    assert(fy->neigh_cells[i]->level >= c->level);
                }

                Face* fz = c->external_faces_r[i];
                if (flag3[i] == false) {
                    if (fz->neigh_cells[i] != c && fz->neigh_cells[i] != NULL) {
                        map<unsigned int, Cell*> temp;
                        temp = refinement(fz->neigh_cells[i]);
                        split_cells.insert(temp.begin(), temp.end());
                    }
                    assert(fz->neigh_cells[i]->level >= c->level);
                }
            }
        }
    } while (flag == false);

    int child_level = c->level + 1;
    int max_level_previous = cells.size() - 1;
    if (child_level > max_level_previous) {
        cells.push_back(vector<Cell*>(0));
        faces.push_back(vector<Face*>(0));
        num_cell.push_back(0);
        num_face.push_back(0);
    }
    // cout << "CD" << endl;
    if (flag) {
        c->isleaf = false;

        double r[3], theta[3], phi[3];
        double dr, dtheta, dphi;
        assert(c->_r[0] < c->_r[1]);
        assert(c->_theta[0] < c->_theta[1]);
        assert(c->_theta[0] < c->_theta[1]);
        r[0] = c->_r[0];
        theta[0] = c->_theta[0];
        phi[0] = c->_phi[0];
        c->get_center(r[1], theta[1], phi[1]);
        c->get_size(dr, dtheta, dphi);
        r[2] = c->_r[1];
        theta[2] = c->_theta[1];
        phi[2] = c->_phi[1];

        Cell* c_child[8];
        int id = 0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    c_child[id] =
                        new Cell(r[k], theta[i], phi[j], r[k + 1], theta[i + 1],
                                 phi[j + 1], child_level, n_parameters, true);
                    assert(c->parameters.size() == n_parameters);
                    assert(c->parameters.size() ==
                           c_child[id]->parameters.size());
                    for (int ip = 0; ip < n_parameters; ip++) {
                        c_child[id]->set_parameter(c->get_parameter(ip), ip);
                    }
                    id++;
                }
            }
        }

        for (int i = 0; i < 8; i++) {
            c->child_cells[i] = c_child[i];
            c_child[i]->set_ordering_forward(c->get_ordering_forward());
            cells[child_level].push_back(c_child[i]);
        }

        Face* internal_child_theta[4];
        Face* internal_child_phi[4];
        Face* internal_child_r[4];

        internal_child_theta[0] =
            new Face(c_child[0], c_child[4], r[1] - 0.25 * dr, theta[1],
                     phi[1] - 0.25 * dphi, NORTH_SOUTH, child_level, true);
        internal_child_theta[1] =
            new Face(c_child[1], c_child[5], r[1] + 0.25 * dr, theta[1],
                     phi[1] - 0.25 * dphi, NORTH_SOUTH, child_level, true);
        internal_child_theta[2] =
            new Face(c_child[2], c_child[6], r[1] - 0.25 * dr, theta[1],
                     phi[1] + 0.25 * dphi, NORTH_SOUTH, child_level, true);
        internal_child_theta[3] =
            new Face(c_child[3], c_child[7], r[1] + 0.25 * dr, theta[1],
                     phi[1] + 0.25 * dphi, NORTH_SOUTH, child_level, true);

        c->set_internal_faces(internal_child_theta, NORTH_SOUTH);

        for (int i = 0; i < 4; i++) {
            faces[child_level].push_back(internal_child_theta[i]);
            leaf_faces.push_back(internal_child_theta[i]);
        }

        internal_child_phi[0] = new Face(
            c_child[0], c_child[2], r[1] - 0.25 * dr, theta[1] - 0.25 * dtheta,
            phi[1], WEST_EAST, child_level, true);
        internal_child_phi[1] = new Face(
            c_child[1], c_child[3], r[1] + 0.25 * dr, theta[1] - 0.25 * dtheta,
            phi[1], WEST_EAST, child_level, true);
        internal_child_phi[2] = new Face(
            c_child[4], c_child[6], r[1] - 0.25 * dr, theta[1] + 0.25 * dtheta,
            phi[1], WEST_EAST, child_level, true);
        internal_child_phi[3] = new Face(
            c_child[5], c_child[7], r[1] + 0.25 * dr, theta[1] + 0.25 * dtheta,
            phi[1], WEST_EAST, child_level, true);

        c->set_internal_faces(internal_child_phi, WEST_EAST);

        for (int i = 0; i < 4; i++) {
            faces[child_level].push_back(internal_child_phi[i]);
            leaf_faces.push_back(internal_child_phi[i]);
        }

        internal_child_r[0] =
            new Face(c_child[0], c_child[1], r[1], theta[1] - 0.25 * dtheta,
                     phi[1] - 0.25 * dphi, RADIUS, child_level, true);
        internal_child_r[1] =
            new Face(c_child[2], c_child[3], r[1], theta[1] - 0.25 * dtheta,
                     phi[1] + 0.25 * dphi, RADIUS, child_level, true);
        internal_child_r[2] =
            new Face(c_child[4], c_child[5], r[1], theta[1] + 0.25 * dtheta,
                     phi[1] - 0.25 * dphi, RADIUS, child_level, true);
        internal_child_r[3] =
            new Face(c_child[6], c_child[7], r[1], theta[1] + 0.25 * dtheta,
                     phi[1] + 0.25 * dphi, RADIUS, child_level, true);

        c->set_internal_faces(internal_child_r, RADIUS);

        for (int i = 0; i < 4; i++) {
            faces[child_level].push_back(internal_child_r[i]);
            leaf_faces.push_back(internal_child_r[i]);
        }
        num_face[child_level] = num_face[child_level] + 12;

        Face* f;
        for (int i = 0; i < 2; i++) {
            f = c->external_faces_theta[i];
            Face* f_child[4];
            if (f->isleaf == true) {
                f->isleaf = false;
                if (i == 1 && f->neigh_cells[1] == NULL) {
                    f_child[0] = new Face(
                        c_child[0 + i * 4], NULL, f->rc - 0.25 * dr, f->thetac,
                        f->phic - 0.25 * dphi, NORTH_SOUTH, child_level, true);
                    f_child[1] = new Face(
                        c_child[1 + i * 4], NULL, f->rc + 0.25 * dr, f->thetac,
                        f->phic - 0.25 * dphi, NORTH_SOUTH, child_level, true);
                    f_child[2] = new Face(
                        c_child[2 + i * 4], NULL, f->rc - 0.25 * dr, f->thetac,
                        f->phic + 0.25 * dphi, NORTH_SOUTH, child_level, true);
                    f_child[3] = new Face(
                        c_child[3 + i * 4], NULL, f->rc + 0.25 * dr, f->thetac,
                        f->phic + 0.25 * dphi, NORTH_SOUTH, child_level, true);
                } else {
                    f_child[0] = new Face(f->neigh_cells[i], c_child[0 + i * 4],
                                          f->rc - 0.25 * dr, f->thetac,
                                          f->phic - 0.25 * dphi, NORTH_SOUTH,
                                          child_level, true);
                    f_child[1] = new Face(f->neigh_cells[i], c_child[1 + i * 4],
                                          f->rc + 0.25 * dr, f->thetac,
                                          f->phic - 0.25 * dphi, NORTH_SOUTH,
                                          child_level, true);
                    f_child[2] = new Face(f->neigh_cells[i], c_child[2 + i * 4],
                                          f->rc - 0.25 * dr, f->thetac,
                                          f->phic + 0.25 * dphi, NORTH_SOUTH,
                                          child_level, true);
                    f_child[3] = new Face(f->neigh_cells[i], c_child[3 + i * 4],
                                          f->rc + 0.25 * dr, f->thetac,
                                          f->phic + 0.25 * dphi, NORTH_SOUTH,
                                          child_level, true);
                }

                f->set_child_faces(f_child[0], f_child[1], f_child[2],
                                   f_child[3]);

                for (int j = 0; j < 4; j++) {
                    faces[child_level].push_back(f_child[j]);
                }
                num_face[child_level] = num_face[child_level] + 4;

                vector<Face*>::iterator it_to_be_deleted =
                    find(leaf_faces.begin(), leaf_faces.end(), f);
                vector<Face*>::iterator it_insert_point =
                    leaf_faces.erase(it_to_be_deleted);
                leaf_faces.insert(it_insert_point, {f_child[0], f_child[1],
                                                    f_child[2], f_child[3]});
            } else {
                f->child_faces[0]->set_neigh_cells(
                    f->child_faces[0]->neigh_cells[i], c_child[0 + i * 4],
                    NORTH_SOUTH);
                f->child_faces[1]->set_neigh_cells(
                    f->child_faces[1]->neigh_cells[i], c_child[1 + i * 4],
                    NORTH_SOUTH);
                f->child_faces[2]->set_neigh_cells(
                    f->child_faces[2]->neigh_cells[i], c_child[2 + i * 4],
                    NORTH_SOUTH);
                f->child_faces[3]->set_neigh_cells(
                    f->child_faces[3]->neigh_cells[i], c_child[3 + i * 4],
                    NORTH_SOUTH);
            }
            c_child[0 + i * 4]->set_external_faces(
                f->child_faces[0], c->internal_faces_theta[0], NORTH_SOUTH);
            c_child[1 + i * 4]->set_external_faces(
                f->child_faces[1], c->internal_faces_theta[1], NORTH_SOUTH);
            c_child[2 + i * 4]->set_external_faces(
                f->child_faces[2], c->internal_faces_theta[2], NORTH_SOUTH);
            c_child[3 + i * 4]->set_external_faces(
                f->child_faces[3], c->internal_faces_theta[3], NORTH_SOUTH);
        }

        for (int i = 0; i < 2; i++) {
            f = c->external_faces_phi[i];
            Face* f_child[4];
            if (f->isleaf) {
                f->isleaf = false;
                if (i == 1 && f->neigh_cells[i] == NULL) {
                    f_child[0] =
                        new Face(c_child[0 + i * 2], NULL, f->rc - 0.25 * dr,
                                 f->thetac - 0.25 * dtheta, f->phic, WEST_EAST,
                                 child_level, true);
                    f_child[1] =
                        new Face(c_child[1 + i * 2], NULL, f->rc + 0.25 * dr,
                                 f->thetac - 0.25 * dtheta, f->phic, WEST_EAST,
                                 child_level, true);
                    f_child[2] =
                        new Face(c_child[4 + i * 2], NULL, f->rc - 0.25 * dr,
                                 f->thetac + 0.25 * dtheta, f->phic, WEST_EAST,
                                 child_level, true);
                    f_child[3] =
                        new Face(c_child[5 + i * 2], NULL, f->rc + 0.25 * dr,
                                 f->thetac + 0.25 * dtheta, f->phic, WEST_EAST,
                                 child_level, true);
                } else {
                    f_child[0] =
                        new Face(f->neigh_cells[i], c_child[0 + i * 2],
                                 f->rc - 0.25 * dr, f->thetac - 0.25 * dtheta,
                                 f->phic, WEST_EAST, child_level, true);
                    f_child[1] =
                        new Face(f->neigh_cells[i], c_child[1 + i * 2],
                                 f->rc + 0.25 * dr, f->thetac - 0.25 * dtheta,
                                 f->phic, WEST_EAST, child_level, true);
                    f_child[2] =
                        new Face(f->neigh_cells[i], c_child[4 + i * 2],
                                 f->rc - 0.25 * dr, f->thetac + 0.25 * dtheta,
                                 f->phic, WEST_EAST, child_level, true);
                    f_child[3] =
                        new Face(f->neigh_cells[i], c_child[5 + i * 2],
                                 f->rc + 0.25 * dr, f->thetac + 0.25 * dtheta,
                                 f->phic, WEST_EAST, child_level, true);
                }

                f->set_child_faces(f_child[0], f_child[1], f_child[2],
                                   f_child[3]);
                for (int j = 0; j < 4; j++) {
                    faces[child_level].push_back(f_child[j]);
                }
                num_face[child_level] = num_face[child_level] + 4;

                vector<Face*>::iterator it_to_be_deleted =
                    find(leaf_faces.begin(), leaf_faces.end(), f);
                vector<Face*>::iterator it_insert_point =
                    leaf_faces.erase(it_to_be_deleted);
                leaf_faces.insert(it_insert_point, {f_child[0], f_child[1],
                                                    f_child[2], f_child[3]});
            } else {
                f->child_faces[0]->set_neigh_cells(
                    f->child_faces[0]->neigh_cells[i], c_child[0 + i * 2],
                    WEST_EAST);
                f->child_faces[1]->set_neigh_cells(
                    f->child_faces[1]->neigh_cells[i], c_child[1 + i * 2],
                    WEST_EAST);
                f->child_faces[2]->set_neigh_cells(
                    f->child_faces[2]->neigh_cells[i], c_child[4 + i * 2],
                    WEST_EAST);
                f->child_faces[3]->set_neigh_cells(
                    f->child_faces[3]->neigh_cells[i], c_child[5 + i * 2],
                    WEST_EAST);
            }
            c_child[0 + i * 2]->set_external_faces(
                f->child_faces[0], c->internal_faces_phi[0], WEST_EAST);
            c_child[1 + i * 2]->set_external_faces(
                f->child_faces[1], c->internal_faces_phi[1], WEST_EAST);
            c_child[4 + i * 2]->set_external_faces(
                f->child_faces[2], c->internal_faces_phi[2], WEST_EAST);
            c_child[5 + i * 2]->set_external_faces(
                f->child_faces[3], c->internal_faces_phi[3], WEST_EAST);
        }

        for (int i = 0; i < 2; i++) {
            f = c->external_faces_r[i];
            Face* f_child[4];
            if (f->isleaf) {
                f->isleaf = false;
                if (i == 1 && f->neigh_cells[i] == NULL) {
                    f_child[0] = new Face(c_child[0 + i * 1], NULL, f->rc,
                                          f->thetac - 0.25 * dtheta,
                                          f->phic - 0.25 * dphi, RADIUS,
                                          child_level, true);
                    f_child[1] = new Face(c_child[2 + i * 1], NULL, f->rc,
                                          f->thetac - 0.25 * dtheta,
                                          f->phic + 0.25 * dphi, RADIUS,
                                          child_level, true);
                    f_child[2] = new Face(c_child[4 + i * 1], NULL, f->rc,
                                          f->thetac + 0.25 * dtheta,
                                          f->phic - 0.25 * dphi, RADIUS,
                                          child_level, true);
                    f_child[3] = new Face(c_child[6 + i * 1], NULL, f->rc,
                                          f->thetac + 0.25 * dtheta,
                                          f->phic + 0.25 * dphi, RADIUS,
                                          child_level, true);
                } else {
                    if (f->neigh_cells[0] != NULL &&
                        f->neigh_cells[1] != NULL) {
                        if (!(abs(f->neigh_cells[0]->_r[1] - f->rc) < 1e-10)) {
                            f->display();
                        }

                        assert(abs(f->neigh_cells[0]->_r[1] -
                                   f->neigh_cells[1]->_r[0]) < 1e-10);
                        assert(abs(f->neigh_cells[0]->_r[1] - f->rc) < 1e-10);
                    }

                    f_child[0] = new Face(f->neigh_cells[i], c_child[0 + i * 1],
                                          f->rc, f->thetac - 0.25 * dtheta,
                                          f->phic - 0.25 * dphi, RADIUS,
                                          child_level, true);
                    f_child[1] = new Face(f->neigh_cells[i], c_child[2 + i * 1],
                                          f->rc, f->thetac - 0.25 * dtheta,
                                          f->phic + 0.25 * dphi, RADIUS,
                                          child_level, true);
                    f_child[2] = new Face(f->neigh_cells[i], c_child[4 + i * 1],
                                          f->rc, f->thetac + 0.25 * dtheta,
                                          f->phic - 0.25 * dphi, RADIUS,
                                          child_level, true);
                    f_child[3] = new Face(f->neigh_cells[i], c_child[6 + i * 1],
                                          f->rc, f->thetac + 0.25 * dtheta,
                                          f->phic + 0.25 * dphi, RADIUS,
                                          child_level, true);
                }

                f->set_child_faces(f_child[0], f_child[1], f_child[2],
                                   f_child[3]);
                for (int j = 0; j < 4; j++) {
                    faces[child_level].push_back(f_child[j]);
                }
                num_face[child_level] = num_face[child_level] + 4;

                vector<Face*>::iterator it_to_be_deleted =
                    find(leaf_faces.begin(), leaf_faces.end(), f);
                vector<Face*>::iterator it_insert_point =
                    leaf_faces.erase(it_to_be_deleted);
                leaf_faces.insert(it_insert_point, {f_child[0], f_child[1],
                                                    f_child[2], f_child[3]});
            } else {
                f->child_faces[0]->set_neigh_cells(
                    f->child_faces[0]->neigh_cells[i], c_child[0 + i * 1],
                    RADIUS);
                f->child_faces[1]->set_neigh_cells(
                    f->child_faces[1]->neigh_cells[i], c_child[2 + i * 1],
                    RADIUS);
                f->child_faces[2]->set_neigh_cells(
                    f->child_faces[2]->neigh_cells[i], c_child[4 + i * 1],
                    RADIUS);
                f->child_faces[3]->set_neigh_cells(
                    f->child_faces[3]->neigh_cells[i], c_child[6 + i * 1],
                    RADIUS);
            }
            c_child[0 + i * 1]->set_external_faces(
                f->child_faces[0], c->internal_faces_r[0], RADIUS);
            c_child[2 + i * 1]->set_external_faces(
                f->child_faces[1], c->internal_faces_r[1], RADIUS);
            c_child[4 + i * 1]->set_external_faces(
                f->child_faces[2], c->internal_faces_r[2], RADIUS);
            c_child[6 + i * 1]->set_external_faces(
                f->child_faces[3], c->internal_faces_r[3], RADIUS);
        }

        vector<Cell*>::iterator c_iter =
            std::find(leaf_cells.begin(), leaf_cells.end(), c);
        // assert(c_iter != leaf_cells.end());
        vector<Cell*>::iterator iter = leaf_cells.erase(c_iter);
        // vector<Cell*>::iterator iter2 =
        leaf_cells.insert(
            iter, {c_child[0], c_child[1], c_child[2], c_child[3], c_child[4],
                   c_child[5], c_child[6], c_child[7]});
        num_leaf_cells = leaf_cells.size();
        num_leaf_faces = leaf_faces.size();

        num_cell[child_level] = num_cell[child_level] + 8;

        vector<Cell*>::iterator c_iter_ro =
            std::find(ro_leaf_cells.begin(), ro_leaf_cells.end(), c);
        vector<Cell*>::iterator iter_ro = ro_leaf_cells.erase(c_iter_ro);
        if (c->get_ordering_forward()) {
            ro_leaf_cells.insert(
                iter_ro, {c_child[0], c_child[2], c_child[6], c_child[4],
                          c_child[5], c_child[7], c_child[3], c_child[1]});
        } else {
            ro_leaf_cells.insert(
                iter_ro, {c_child[1], c_child[3], c_child[7], c_child[5],
                          c_child[4], c_child[6], c_child[2], c_child[0]});
        }

        // assert(iter2 != leaf_cells.end());
        // split_cells[c->id]=iter2;
        split_cells.insert(pair<unsigned int, Cell*>(c->id, c));
        // cout << "id=" << c->id << ", " << (*split_cells[c->id])->isleaf <<
        // endl; return iter2;
    }

    return split_cells;
    // else
    // {
    //     cout << "leaf_cells.end()" << endl;
    //     return leaf_cells.end();
    // }
}

map<unsigned int, Cell*> Mesh::refinement(int index_to_be_refined) {
    Cell* c = leaf_cells[index_to_be_refined];
    map<unsigned int, Cell*> map_cells = this->refinement(c);
    return map_cells;
}

void Mesh::get_model_parameter_from_mesh(VectorXd& m, int ith) {
    m.resize(leaf_cells.size());
#pragma omp parallel for
    for (int i = 0; i < leaf_cells.size(); i++) {
        m(i) = leaf_cells[i]->get_parameter(ith);
    }
}

void Mesh::rearrange_id() {
#pragma omp parallel for
    for (int i = 0; i < leaf_cells.size(); i++) {
        leaf_cells[i]->set_id(i);
    }
    new_to_old_index.resize(leaf_cells.size());
    old_to_new_index.resize(leaf_cells.size());
#pragma omp parallel for
    for (int i = 0; i < leaf_cells.size(); i++) {
        new_to_old_index[i] = ro_leaf_cells[i]->get_id();
    }
#pragma omp parallel for
    for (int i = 0; i < leaf_cells.size(); i++) {
        int old_index = new_to_old_index[i];
        int new_index = i;
        old_to_new_index[old_index] = new_index;
    }
}

bool Mesh::great_equal(long double left, long double right) {
    bool temp = false;
    // left>=right?
    long double a = left - right;
    if (a > 0.) temp = true;             // a>0, OK.
    if (std::abs(a) < TOL) temp = true;  // a=0, OK.
    return temp;
}

void Mesh::rthetaphi_xyz(long double r, long double theta, long double phi,
                         long double& x, long double& y, long double& z) {
    // we must check the range of r, theta and phi
    assert(great_equal(r, 0.));
    assert(great_equal(theta, 0.));
    assert(great_equal(GS::PI, theta));
    assert(great_equal(phi, 0.));
    // cout<<phi<<endl;
    assert(great_equal(2.0 * GS::PI, phi));

    x = std::sin(theta) * std::cos(phi) * r;
    y = std::sin(theta) * std::sin(phi) * r;
    z = std::cos(theta) * r;
    return;
}
void Mesh::sort(int level) {
    for (int i = 0; i < cells[level].size(); i++) {
        Cell* c = cells[level][i];
        // sort faces
        if (c->external_faces_r[0]->rc > c->external_faces_r[1]->rc) {
            Face* t = c->external_faces_r[0];
            c->external_faces_r[0] = c->external_faces_r[1];
            c->external_faces_r[1] = t;
        }

        if (c->external_faces_theta[0]->thetac >
            c->external_faces_theta[1]->thetac) {
            Face* t = c->external_faces_theta[0];
            c->external_faces_theta[0] = c->external_faces_theta[1];
            c->external_faces_theta[1] = t;
        }

        if (c->external_faces_phi[0]->phic > c->external_faces_phi[1]->phic) {
            Face* t = c->external_faces_phi[0];
            c->external_faces_phi[0] = c->external_faces_phi[1];
            c->external_faces_phi[1] = t;
        }

        // sort neighbouring cells
        for (int j = 0; j < 2; j++) {
            Cell *n1, *n2;
            n1 = c->external_faces_r[j]->neigh_cells[0];
            n2 = c->external_faces_r[j]->neigh_cells[1];

            double rc1, thetac1, phic1;
            double rc2, thetac2, phic2;

            if (n1 != NULL && n2 != NULL) {
                n1->get_center(rc1, thetac1, phic1);
                n2->get_center(rc2, thetac2, phic2);
                if (rc1 > rc2) {
                    c->external_faces_r[j]->neigh_cells[0] = n2;
                    c->external_faces_r[j]->neigh_cells[1] = n1;
                }
            } else {
                double rc0, thetac0, phic0;
                c->get_center(rc0, thetac0, phic0);
                if (rc0 > (c->external_faces_r[j]->rc)) {
                    if (n2 == NULL && n1 == c) {
                        c->external_faces_r[j]->neigh_cells[0] = NULL;
                        c->external_faces_r[j]->neigh_cells[1] = c;
                    }
                } else {
                    if (n1 == NULL && n2 == c) {
                        c->external_faces_r[j]->neigh_cells[0] = c;
                        c->external_faces_r[j]->neigh_cells[1] = NULL;
                    }
                }
            }
        }

        for (int j = 0; j < 2; j++) {
            Cell *n1, *n2;
            n1 = c->external_faces_theta[j]->neigh_cells[0];
            n2 = c->external_faces_theta[j]->neigh_cells[1];

            double rc1, thetac1, phic1;
            double rc2, thetac2, phic2;

            if (n1 != NULL && n2 != NULL) {
                n1->get_center(rc1, thetac1, phic1);
                n2->get_center(rc2, thetac2, phic2);
                if (thetac1 > thetac2) {
                    c->external_faces_theta[j]->neigh_cells[0] = n2;
                    c->external_faces_theta[j]->neigh_cells[1] = n1;
                }
            } else {
                double rc0, thetac0, phic0;
                c->get_center(rc0, thetac0, phic0);
                if (thetac0 > (c->external_faces_theta[j]->thetac)) {
                    if (n2 == NULL && n1 == c) {
                        c->external_faces_theta[j]->neigh_cells[0] = NULL;
                        c->external_faces_theta[j]->neigh_cells[1] = c;
                    }
                } else {
                    if (n1 == NULL && n2 == c) {
                        c->external_faces_theta[j]->neigh_cells[0] = c;
                        c->external_faces_theta[j]->neigh_cells[1] = NULL;
                    }
                }
            }
        }

        for (int j = 0; j < 2; j++) {
            Cell *n1, *n2;
            n1 = c->external_faces_phi[j]->neigh_cells[0];
            n2 = c->external_faces_phi[j]->neigh_cells[1];

            double rc1, thetac1, phic1;
            double rc2, thetac2, phic2;

            if (n1 != NULL && n2 != NULL) {
                n1->get_center(rc1, thetac1, phic1);
                n2->get_center(rc2, thetac2, phic2);
                if (phic1 > phic2) {
                    c->external_faces_phi[j]->neigh_cells[0] = n2;
                    c->external_faces_phi[j]->neigh_cells[1] = n1;
                }
            } else {
                double rc0, thetac0, phic0;
                c->get_center(rc0, thetac0, phic0);
                if (phic0 > (c->external_faces_phi[j]->phic)) {
                    if (n2 == NULL && n1 == c) {
                        c->external_faces_phi[j]->neigh_cells[0] = NULL;
                        c->external_faces_phi[j]->neigh_cells[1] = c;
                    }
                } else {
                    if (n1 == NULL && n2 == c) {
                        c->external_faces_phi[j]->neigh_cells[0] = c;
                        c->external_faces_phi[j]->neigh_cells[1] = NULL;
                    }
                }
            }
        }
    }
}

void Mesh::fill_data(int offset_i, int offset_j, int offset_k, double*** data,
                     Cell* c, int max_level, int ith_para) {
    if (c->isleaf) {
        int level_difference = max_level - c->level;
        int N = pow(2, level_difference);
        double value = c->get_parameter(ith_para);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    data[offset_i + i][offset_j + j][offset_k + k] = value;
                }
            }
        }
    } else {
        int level_difference = max_level - c->child_cells[0]->level;
        int N = pow(2, level_difference);
        //递归
        fill_data(offset_i, offset_j, offset_k, data, c->child_cells[0],
                  max_level, ith_para);
        fill_data(offset_i, offset_j, offset_k + N, data, c->child_cells[1],
                  max_level, ith_para);
        fill_data(offset_i, offset_j + N, offset_k, data, c->child_cells[2],
                  max_level, ith_para);
        fill_data(offset_i, offset_j + N, offset_k + N, data, c->child_cells[3],
                  max_level, ith_para);
        fill_data(offset_i + N, offset_j, offset_k, data, c->child_cells[4],
                  max_level, ith_para);
        fill_data(offset_i + N, offset_j, offset_k + N, data, c->child_cells[5],
                  max_level, ith_para);
        fill_data(offset_i + N, offset_j + N, offset_k, data, c->child_cells[6],
                  max_level, ith_para);
        fill_data(offset_i + N, offset_j + N, offset_k + N, data,
                  c->child_cells[7], max_level, ith_para);
    }
}

int Mesh::out_model_netcdf(string filename, int ith_para, string VAL_NAME,
                           string VAL_UNITS) {
    int max_level = cells.size() - 1;
    int N = std::pow(2, max_level);
    int NR = nr * N;
    int NLAT = ntheta * N;
    int NLON = nphi * N;

    assert(cells[0].size() == nr * ntheta * nphi);

    int n_total = std::pow(8, max_level) * nr * ntheta * nphi;
    double*** DENSITY_DATA;
    DENSITY_DATA = new double**[NLAT];
    for (unsigned int i = 0; i < NLAT; i++) {
        DENSITY_DATA[i] = new double*[NLON];
        for (unsigned int j = 0; j < NLON; j++) {
            DENSITY_DATA[i][j] = new double[NR];
        }
    }

    for (int i = 0; i < ntheta; i++) {
        for (int j = 0; j < nphi; j++) {
            for (int k = 0; k < nr; k++) {
                fill_data(i * N, j * N, k * N, DENSITY_DATA,
                          cells[0][i * nphi * nr + j * nr + k], max_level,
                          ith_para);
            }
        }
    }

    // writh data to netcdf file
    double* lats = new double[NLAT];
    double* lons = new double[NLON];
    double* rs = new double[NR];

    double** lats_bnd = new double*[NLAT];
    double** lons_bnd = new double*[NLON];
    double** r_bnd = new double*[NR];

    // lats=new double[_lats_cells];
    // lons=new double[_long_cells];
    // rs=new double [NR0];

    double theta_space = (theta_lim[1] - theta_lim[0]) / NLAT;
    for (unsigned int i = 0; i < NLAT; i++) {
        lats[i] = theta_lim[0] + 0.5 * theta_space + i * theta_space;
        lats[i] = 90.0 - lats[i] * 180.0 / GS::PI;
        lats_bnd[i] = new double[2];
        lats_bnd[i][0] =
            90.0 - (theta_lim[0] + i * theta_space) * 180.0 / GS::PI;
        lats_bnd[i][1] =
            90.0 - (theta_lim[0] + (i + 1) * theta_space) * 180.0 / GS::PI;
    }

    double phi_space = (phi_lim[1] - phi_lim[0]) / NLON;
    for (unsigned int j = 0; j < NLON; j++) {
        lons[j] = phi_lim[0] + 0.5 * phi_space + j * phi_space;
        lons[j] = lons[j] * 180.0 / GS::PI - 180.0;
        lons_bnd[j] = new double[2];

        lons_bnd[j][0] = (phi_lim[0] + j * phi_space) * 180.0 / GS::PI - 180.0;
        lons_bnd[j][1] =
            (phi_lim[0] + (j + 1) * phi_space) * 180.0 / GS::PI - 180.0;
    }

    // double r_space = (r_lim[1] - r_lim[0]) / NR;
    // for (unsigned int k = 0; k < NR; k++)
    // {
    //     rs[k] = r_lim[0] + 0.5 * r_space + k * r_space;
    //     r_bnd[k] = new double[2];
    //     r_bnd[k][0] = r_lim[0] + k * r_space;
    //     r_bnd[k][1] = r_lim[0] + (k + 1) * r_space;
    // }
    for (unsigned int k = 0; k < nr; k++) {
        double r_space = (r_points(k + 1) - r_points(k)) / (1.0 * N);
        for (unsigned int k2 = 0; k2 < N; k2++) {
            int index = k * N + k2;
            rs[index] = r_points(k) + 0.5 * r_space + k2 * r_space;
            r_bnd[index] = new double[2];
            r_bnd[index][0] = r_points(k) + k2 * r_space;
            r_bnd[index][1] = r_points(k) + (k2 + 1) * r_space;
        }
    }

    // Names of things.
    const char* LAT_NAME = "latitude";
    const char* LON_NAME = "longitude";
    const char* VAL_NAME_C = VAL_NAME.c_str();
    const char* R_NAME = "radius";

    string UNITS = "units";
    string POSITIVE = "positive";
    string DEGREES_EAST = "degrees_east";
    string DEGREES_NORTH = "degrees_north";
    string METER = "meter";
    string UP = "up";
    // For the units attributes.

    string LAT_UNITS = "degrees_north";
    string LON_UNITS = "degrees_east";
    string BOUNDS = "bounds";
    string LATBND_NAME = "latbnd";
    string LONBND_NAME = "lonbnd";
    string RBND_NAME = "rbnd";

    try {
        // Create the file. The Replace parameter tells netCDF to overwrite
        // this file, if it already exists.
        NcFile test(filename, NcFile::replace);
        test.putAtt("Conventions", "CF-1.7");
        int pixel_registration[1] = {1};
        string NODE_OFFSET = "node_offset";
        // test.putAtt(NODE_OFFSET, "1");
        // test.putAtt(NODE_OFFSET, ncInt, 1, pixel_registration);
        // Define the dimensions. NetCDF will hand back an ncDim object for
        // each.
        NcDim latDim = test.addDim(LAT_NAME, NLAT);
        NcDim lonDim = test.addDim(LON_NAME, NLON);
        NcDim rDim = test.addDim(R_NAME, NR);

        NcVar latVar = test.addVar(LAT_NAME, ncDouble, latDim);
        NcVar lonVar = test.addVar(LON_NAME, ncDouble, lonDim);
        NcVar rVar = test.addVar(R_NAME, ncDouble, rDim);

        // Define units attributes for coordinate vars. This attaches a
        // text attribute to each of the coordinate variables, containing
        // the units.
        latVar.putAtt("long_name", "Latitude");
        latVar.putAtt("standard_name", "Latitude");
        latVar.putAtt(UNITS, DEGREES_NORTH);
        lonVar.putAtt(UNITS, DEGREES_EAST);
        lonVar.putAtt("long_name", "Longitude");
        lonVar.putAtt("standard_name", "Longitude");
        rVar.putAtt(UNITS, METER);
        rVar.putAtt(POSITIVE, UP);

        // Write the coordinate variable data to the file.
        latVar.putVar(lats);
        lonVar.putVar(lons);
        rVar.putVar(rs);

        // bounds for cells
        NcDim bndDim = test.addDim("bnd", 2);
        vector<NcDim> dimVector_latbnd;
        dimVector_latbnd.push_back(latDim);
        dimVector_latbnd.push_back(bndDim);
        NcVar latbnd = test.addVar(LATBND_NAME, ncDouble, dimVector_latbnd);
        latVar.putAtt(BOUNDS, LATBND_NAME);

        vector<NcDim> dimVector_lonbnd;
        dimVector_lonbnd.push_back(lonDim);
        dimVector_lonbnd.push_back(bndDim);
        NcVar lonbnd = test.addVar(LONBND_NAME, ncDouble, dimVector_lonbnd);
        lonVar.putAtt(BOUNDS, LONBND_NAME);

        vector<NcDim> dimVector_rbnd;
        dimVector_rbnd.push_back(rDim);
        dimVector_rbnd.push_back(bndDim);
        NcVar rbnd = test.addVar(RBND_NAME, ncDouble, dimVector_rbnd);
        rVar.putAtt(BOUNDS, RBND_NAME);

        vector<size_t> startp_bnd, countp_bnd;
        startp_bnd.push_back(0);
        startp_bnd.push_back(0);
        countp_bnd.push_back(1);
        countp_bnd.push_back(1);
        for (int i = 0; i < NLAT; i++) {
            startp_bnd[0] = i;
            for (int j = 0; j < 2; j++) {
                startp_bnd[1] = j;
                double a = lats_bnd[i][j];
                latbnd.putVar(startp_bnd, countp_bnd, &a);
            }
        }

        for (int i = 0; i < NLON; i++) {
            startp_bnd[0] = i;
            for (int j = 0; j < 2; j++) {
                startp_bnd[1] = j;
                double a = lons_bnd[i][j];
                lonbnd.putVar(startp_bnd, countp_bnd, &a);
            }
        }

        for (int i = 0; i < NR; i++) {
            startp_bnd[0] = i;
            for (int j = 0; j < 2; j++) {
                startp_bnd[1] = j;
                double a = r_bnd[i][j];
                rbnd.putVar(startp_bnd, countp_bnd, &a);
            }
        }

        // Define the netCDF variables for the density data
        vector<NcDim> dimVector;
        dimVector.push_back(latDim);
        dimVector.push_back(lonDim);
        dimVector.push_back(rDim);
        NcVar denVar = test.addVar(VAL_NAME_C, ncDouble, dimVector);

        denVar.putAtt(UNITS, VAL_UNITS);
        // denVar.putAtt(NODE_OFFSET, "1");

        // Write the density data;
        vector<size_t> startp, countp;
        startp.push_back(0);
        startp.push_back(0);
        startp.push_back(0);

        countp.push_back(1);
        countp.push_back(1);
        countp.push_back(1);

        for (size_t i = 0; i < NLAT; i++) {
            for (size_t j = 0; j < NLON; j++) {
                for (size_t k = 0; k < NR; k++) {
                    startp[0] = i;
                    startp[1] = j;
                    startp[2] = k;
                    double a = DENSITY_DATA[i][j][k];
                    denVar.putVar(startp, countp, &a);
                }
            }
        }

        // free resources
        for (int i = 0; i < NLAT; i++) {
            for (int j = 0; j < NLON; j++) {
                delete[] DENSITY_DATA[i][j];
                DENSITY_DATA[i][j] = NULL;
            }
            delete[] DENSITY_DATA[i];
            DENSITY_DATA[i] = NULL;
        }
        delete[] DENSITY_DATA;
        DENSITY_DATA = NULL;

        delete[] lats;
        lats = NULL;

        delete[] lons;
        lons = NULL;

        delete[] rs;
        rs = NULL;

        for (int i = 0; i < NLAT; i++) {
            delete[] lats_bnd[i];
            lats_bnd[i] = NULL;
        }
        delete[] lats_bnd;
        lats_bnd = NULL;
        for (int j = 0; j < NLON; j++) {
            delete[] lons_bnd[j];
            lons_bnd[j] = NULL;
        }
        delete[] lons_bnd;
        lons_bnd = NULL;
        for (int k = 0; k < NR; k++) {
            delete[] r_bnd[k];
            r_bnd[k] = NULL;
        }
        delete[] r_bnd;
        r_bnd = NULL;

        // denVar.putVar(DENSITY_DATA);
        cout << "The model has been written to NetCDF file: " << filename
             << endl;
        return 0;
    } catch (NcException& e) {
        e.what();
        return NC_ERR;
    }
}

void Mesh::refine_cell_across_a_interface(Cell* c, double r_interface,
                                          double r_resolution,
                                          double upper_part, double lower_part,
                                          int i_para) {
    double dep[2] = {c->_r[0], c->_r[1]};
    if (c->isleaf) {
        if (dep[0] < r_interface && dep[1] > r_interface) {
            // cout << "true" << endl;
            // cout << setprecision(10) << dep[0] << ", " << setprecision(10) <<
            // r_interface << ", " << setprecision(10) << dep[1] << "\t\t\t" <<
            // 0.5 * (dep[0] + dep[1]) << "\t" << std::abs(r_interface - 0.5 *
            // (dep[0] + dep[1])) << endl;
            if (r_interface < 0.5 * (dep[0] + dep[1])) {
                c->set_parameter(upper_part, i_para);
            } else {
                c->set_parameter(lower_part, i_para);
            }

            map<unsigned, Cell*> iters = this->refinement(c);
            if (iters.size() != 0) {
                if (std::abs(r_interface - 0.5 * (dep[0] + dep[1])) <
                    r_resolution) {
                    c->child_cells[0]->set_parameter(lower_part, i_para);
                    c->child_cells[2]->set_parameter(lower_part, i_para);
                    c->child_cells[4]->set_parameter(lower_part, i_para);
                    c->child_cells[6]->set_parameter(lower_part, i_para);

                    c->child_cells[1]->set_parameter(upper_part, i_para);
                    c->child_cells[3]->set_parameter(upper_part, i_para);
                    c->child_cells[5]->set_parameter(upper_part, i_para);
                    c->child_cells[7]->set_parameter(upper_part, i_para);
                } else if (r_interface > 0.5 * (dep[0] + dep[1])) {
                    // cout << "1" << endl;
                    assert(c->child_cells[0]->isleaf);
                    assert(c->child_cells[2]->isleaf);
                    assert(c->child_cells[4]->isleaf);
                    assert(c->child_cells[6]->isleaf);
                    // cout << "11" << endl;
                    c->child_cells[0]->set_parameter(lower_part, i_para);
                    c->child_cells[2]->set_parameter(lower_part, i_para);
                    c->child_cells[4]->set_parameter(lower_part, i_para);
                    c->child_cells[6]->set_parameter(lower_part, i_para);

                    // cout << "1" << endl;
                    refine_cell_across_a_interface(
                        c->child_cells[1], r_interface, r_resolution,
                        upper_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[3], r_interface, r_resolution,
                        upper_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[5], r_interface, r_resolution,
                        upper_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[7], r_interface, r_resolution,
                        upper_part, lower_part, i_para);
                } else {
                    c->child_cells[1]->set_parameter(upper_part, i_para);
                    c->child_cells[3]->set_parameter(upper_part, i_para);
                    c->child_cells[5]->set_parameter(upper_part, i_para);
                    c->child_cells[7]->set_parameter(upper_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[0], r_interface, r_resolution,
                        upper_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[2], r_interface, r_resolution,
                        upper_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[4], r_interface, r_resolution,
                        upper_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[6], r_interface, r_resolution,
                        upper_part, lower_part, i_para);
                }
                // cout << "B" << endl;
            }
        }
    } else {
        // cout << "XX" << endl;
        if (std::abs(r_interface - 0.5 * (dep[0] + dep[1])) < r_resolution) {
            c->child_cells[0]->set_parameter(lower_part, i_para);
            c->child_cells[2]->set_parameter(lower_part, i_para);
            c->child_cells[4]->set_parameter(lower_part, i_para);
            c->child_cells[6]->set_parameter(lower_part, i_para);

            c->child_cells[1]->set_parameter(upper_part, i_para);
            c->child_cells[3]->set_parameter(upper_part, i_para);
            c->child_cells[5]->set_parameter(upper_part, i_para);
            c->child_cells[7]->set_parameter(upper_part, i_para);
        } else if (r_interface > 0.5 * (dep[0] + dep[1])) {
            // cout << "1" << endl;
            // cout << "11" << endl;
            c->child_cells[0]->set_parameter(lower_part, i_para);
            c->child_cells[2]->set_parameter(lower_part, i_para);
            c->child_cells[4]->set_parameter(lower_part, i_para);
            c->child_cells[6]->set_parameter(lower_part, i_para);

            // cout << "1" << endl;
            refine_cell_across_a_interface(c->child_cells[1], r_interface,
                                           r_resolution, upper_part, lower_part,
                                           i_para);
            refine_cell_across_a_interface(c->child_cells[3], r_interface,
                                           r_resolution, upper_part, lower_part,
                                           i_para);
            refine_cell_across_a_interface(c->child_cells[5], r_interface,
                                           r_resolution, upper_part, lower_part,
                                           i_para);
            refine_cell_across_a_interface(c->child_cells[7], r_interface,
                                           r_resolution, upper_part, lower_part,
                                           i_para);
        } else {
            c->child_cells[1]->set_parameter(upper_part, i_para);
            c->child_cells[3]->set_parameter(upper_part, i_para);
            c->child_cells[5]->set_parameter(upper_part, i_para);
            c->child_cells[7]->set_parameter(upper_part, i_para);
            refine_cell_across_a_interface(c->child_cells[0], r_interface,
                                           r_resolution, upper_part, lower_part,
                                           i_para);
            refine_cell_across_a_interface(c->child_cells[2], r_interface,
                                           r_resolution, upper_part, lower_part,
                                           i_para);
            refine_cell_across_a_interface(c->child_cells[4], r_interface,
                                           r_resolution, upper_part, lower_part,
                                           i_para);
            refine_cell_across_a_interface(c->child_cells[6], r_interface,
                                           r_resolution, upper_part, lower_part,
                                           i_para);
        }
    }
}

void Mesh::refine_cell_across_2_interfaces(
    Cell* c, double r_interface1, double r_interface2, double r_resolution,
    double upper_part, double middle_part, double lower_part, int i_para) {
    double dep[2] = {c->_r[0], c->_r[1]};
    if (c->isleaf) {
        if (dep[0] < r_interface1 && dep[1] > r_interface2) {
            cout << "true" << endl;
            // cout << setprecision(10) << dep[0] << ", " << setprecision(10) <<
            // r_interface << ", " << setprecision(10) << dep[1] << "\t\t\t" <<
            // 0.5 * (dep[0] + dep[1]) << "\t" << std::abs(r_interface - 0.5 *
            // (dep[0] + dep[1])) << endl;
            if (r_interface2 < 0.5 * (dep[0] + dep[1])) {
                c->set_parameter(upper_part, i_para);
            } else if (r_interface1 > 0.5 * (dep[0] + dep[1])) {
                c->set_parameter(lower_part, i_para);
            } else {
                c->set_parameter(middle_part, i_para);
            }

            map<unsigned, Cell*> iters = this->refinement(c);
            if (iters.size() != 0) {
                if (std::abs(r_interface1 - 0.5 * (dep[0] + dep[1])) <
                    r_resolution) {
                    if (r_interface2 - r_interface1 > r_resolution) {
                        c->child_cells[0]->set_parameter(lower_part, i_para);
                        c->child_cells[2]->set_parameter(lower_part, i_para);
                        c->child_cells[4]->set_parameter(lower_part, i_para);
                        c->child_cells[6]->set_parameter(lower_part, i_para);

                        refine_cell_across_a_interface(
                            c->child_cells[1], r_interface2, r_resolution,
                            upper_part, middle_part, i_para);
                        refine_cell_across_a_interface(
                            c->child_cells[3], r_interface2, r_resolution,
                            upper_part, middle_part, i_para);
                        refine_cell_across_a_interface(
                            c->child_cells[5], r_interface2, r_resolution,
                            upper_part, middle_part, i_para);
                        refine_cell_across_a_interface(
                            c->child_cells[7], r_interface2, r_resolution,
                            upper_part, middle_part, i_para);
                    } else {
                        c->child_cells[0]->set_parameter(lower_part, i_para);
                        c->child_cells[2]->set_parameter(lower_part, i_para);
                        c->child_cells[4]->set_parameter(lower_part, i_para);
                        c->child_cells[6]->set_parameter(lower_part, i_para);

                        c->child_cells[1]->set_parameter(upper_part, i_para);
                        c->child_cells[3]->set_parameter(upper_part, i_para);
                        c->child_cells[5]->set_parameter(upper_part, i_para);
                        c->child_cells[7]->set_parameter(upper_part, i_para);
                    }
                } else if (std::abs(r_interface2 - 0.5 * (dep[0] + dep[1])) <
                           0) {
                    if (r_interface2 - r_interface1 > r_resolution) {
                        c->child_cells[1]->set_parameter(upper_part, i_para);
                        c->child_cells[3]->set_parameter(upper_part, i_para);
                        c->child_cells[5]->set_parameter(upper_part, i_para);
                        c->child_cells[7]->set_parameter(upper_part, i_para);

                        // cout << "1" << endl;
                        refine_cell_across_a_interface(
                            c->child_cells[0], r_interface1, r_resolution,
                            middle_part, lower_part, i_para);
                        refine_cell_across_a_interface(
                            c->child_cells[2], r_interface1, r_resolution,
                            middle_part, lower_part, i_para);
                        refine_cell_across_a_interface(
                            c->child_cells[4], r_interface1, r_resolution,
                            middle_part, lower_part, i_para);
                        refine_cell_across_a_interface(
                            c->child_cells[6], r_interface1, r_resolution,
                            middle_part, lower_part, i_para);
                    } else {
                        c->child_cells[0]->set_parameter(lower_part, i_para);
                        c->child_cells[2]->set_parameter(lower_part, i_para);
                        c->child_cells[4]->set_parameter(lower_part, i_para);
                        c->child_cells[6]->set_parameter(lower_part, i_para);

                        c->child_cells[1]->set_parameter(upper_part, i_para);
                        c->child_cells[3]->set_parameter(upper_part, i_para);
                        c->child_cells[5]->set_parameter(upper_part, i_para);
                        c->child_cells[7]->set_parameter(upper_part, i_para);
                    }
                } else if (r_interface1 > 0.5 * (dep[0] + dep[1])) {
                    c->child_cells[0]->set_parameter(lower_part, i_para);
                    c->child_cells[2]->set_parameter(lower_part, i_para);
                    c->child_cells[4]->set_parameter(lower_part, i_para);
                    c->child_cells[6]->set_parameter(lower_part, i_para);

                    refine_cell_across_2_interfaces(
                        c->child_cells[1], r_interface1, r_interface2,
                        r_resolution, upper_part, middle_part, lower_part,
                        i_para);
                    refine_cell_across_2_interfaces(
                        c->child_cells[3], r_interface1, r_interface2,
                        r_resolution, upper_part, middle_part, lower_part,
                        i_para);
                    refine_cell_across_2_interfaces(
                        c->child_cells[5], r_interface1, r_interface2,
                        r_resolution, upper_part, middle_part, lower_part,
                        i_para);
                    refine_cell_across_2_interfaces(
                        c->child_cells[7], r_interface1, r_interface2,
                        r_resolution, upper_part, middle_part, lower_part,
                        i_para);
                } else if (r_interface2 < 0.5 * (dep[0] + dep[1])) {
                    c->child_cells[1]->set_parameter(upper_part, i_para);
                    c->child_cells[3]->set_parameter(upper_part, i_para);
                    c->child_cells[5]->set_parameter(upper_part, i_para);
                    c->child_cells[7]->set_parameter(upper_part, i_para);

                    refine_cell_across_2_interfaces(
                        c->child_cells[0], r_interface1, r_interface2,
                        r_resolution, upper_part, middle_part, lower_part,
                        i_para);
                    refine_cell_across_2_interfaces(
                        c->child_cells[2], r_interface1, r_interface2,
                        r_resolution, upper_part, middle_part, lower_part,
                        i_para);
                    refine_cell_across_2_interfaces(
                        c->child_cells[4], r_interface1, r_interface2,
                        r_resolution, upper_part, middle_part, lower_part,
                        i_para);
                    refine_cell_across_2_interfaces(
                        c->child_cells[6], r_interface1, r_interface2,
                        r_resolution, upper_part, middle_part, lower_part,
                        i_para);
                } else if (r_interface1 < 0.5 * (dep[0] + dep[1]) &&
                           r_interface2 > 0.5 * (dep[0] + dep[1])) {
                    c->child_cells[1]->set_parameter(upper_part, i_para);
                    c->child_cells[3]->set_parameter(upper_part, i_para);
                    c->child_cells[5]->set_parameter(upper_part, i_para);
                    c->child_cells[7]->set_parameter(upper_part, i_para);

                    refine_cell_across_a_interface(
                        c->child_cells[0], r_interface1, r_resolution,
                        middle_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[2], r_interface1, r_resolution,
                        middle_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[4], r_interface1, r_resolution,
                        middle_part, lower_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[6], r_interface1, r_resolution,
                        middle_part, lower_part, i_para);

                    refine_cell_across_a_interface(
                        c->child_cells[1], r_interface2, r_resolution,
                        upper_part, middle_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[3], r_interface2, r_resolution,
                        upper_part, middle_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[5], r_interface2, r_resolution,
                        upper_part, middle_part, i_para);
                    refine_cell_across_a_interface(
                        c->child_cells[7], r_interface2, r_resolution,
                        upper_part, middle_part, i_para);
                } else {
                    cout << "r_interface1: " << r_interface1 << endl;
                    cout << "r_interface2: " << r_interface2 << endl;
                    cout << "middle depth: " << 0.5 * (dep[0] + dep[1]) << endl;
                    cout << "resolution: " << r_resolution << endl;
                    cout << "B" << endl;
                }
                // cout << "B" << endl;
            }
        } else {
            cout << "C" << endl;
        }
    }
}

void Mesh::show_ordering() {
    int max_level = cells.size() - 1;
    if (max_level == 0) {
        for (int i = 0; i < ntheta; i++) {
            for (int k = 0; k < nr; k++) {
                for (int j = 0; j < nphi; j++) {
                    int id = i * nr * nphi + j * nr + k;
                    int original_id = leaf_cells[id]->get_id();
                    assert(id == original_id);
                    int reordered_id = old_to_new_index[id];
                    double r, theta, phi;
                    leaf_cells[id]->get_center(r, theta, phi);
                    cout << "{" << original_id << ", " << reordered_id
                         << ((leaf_cells[id]->get_ordering_forward())
                                 ? " forward"
                                 : " inverse")
                         << ": (" << r / 1000.0 << ","
                         << 90.0 - theta / GS::PI * 180 << ", "
                         << phi / GS::PI * 180.0 - 180.0 << ")}";
                    cout << "\t";
                }
                cout << endl;
            }
            cout << endl << endl;
        }
    } else {
        cout << "{original_id, new_id: (r, lat, lon)}" << endl;
        for (int i = 0; i < ro_leaf_cells.size(); ++i) {
            int reordered_id = i;
            int original_id = ro_leaf_cells[i]->get_id();
            assert(old_to_new_index[original_id] == i);
            double r, theta, phi;
            ro_leaf_cells[i]->get_center(r, theta, phi);
            cout << "{" << original_id << ", " << reordered_id
                 << ((ro_leaf_cells[i]->get_ordering_forward()) ? (" forward")
                                                                : (" inverse"))
                 << ": (" << r / 1000.0 << "," << 90.0 - theta / GS::PI * 180
                 << ", " << phi / GS::PI * 180.0 - 180.0 << ")}" << endl;
        }
    }
}
