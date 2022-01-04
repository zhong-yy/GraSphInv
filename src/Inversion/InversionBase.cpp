#include "InversionBase.h"
InversionBase::InversionBase()
    : Fwd(),
      a_s(1e-3),
      a_r(1e0),
      a_theta(1e0),
      a_phi(1e0),
      a_crg(1e0),
      target_misfit(1),
      max_lambda(1e6),
      depth_weighting_factor(2),
      n_lambda(20),
      lambda_decreasing_rate(0.5),
      reference_surface(6378137),
      interpolator_m0(NULL),
      interpolator_m0s(NULL),
      use_cross_gradient_constraint(false),
      use_petrophysical_constraint(false) {
    mesh.set_n_parameter(5);
}

InversionBase::InversionBase(const Mesh& mesh_, const Observation& ob_,
                             unsigned long long field_flag_)
    : Fwd(mesh_, ob_, field_flag_),
      target_misfit(1),
      max_lambda(1e6),
      depth_weighting_factor(2),
      n_lambda(20),
      lambda_decreasing_rate(0.5),
      reference_surface(6378137),
      interpolator_m0(NULL),
      interpolator_m0s(NULL),
      use_cross_gradient_constraint(false),
      use_petrophysical_constraint(false) {
    int n_field = field_flag.count();
    this->init_matrices();

    // initialise the reference model as zero
    m0 = VectorXd::Constant(Nm, 1, 0);
    m_min = VectorXd::Constant(Nm, 1, -1e6);
    m_max = VectorXd::Constant(Nm, 1, 1e6);
    m0_s = VectorXd::Constant(Nm, 1, 0);
    // cout<<"A"<<endl;
    Wd.resize(Nd, Nd);
    Wd.setIdentity();

    // cout << "here" << endl;
    m.resize(Nm);

    m_ini.resize(Nm);
    m_ini = VectorXd::Constant(Nm, 1, 0);

    mesh.set_n_parameter(5);

    this->set_depth_weighting(this->depth_weighting_factor);

    constraint_depth[0] = -1e8;
    constraint_depth[1] = 1e8;
    constraint_lat[0] = -1e8;
    constraint_lat[1] = 1e8;
    constraint_lon[0] = -1e8;
    constraint_lon[1] = 1e8;
    // this->update_S_crg();
}

InversionBase::~InversionBase() {
    if (this->interpolator_m0 != NULL) {
        delete interpolator_m0;
    }
    if (this->interpolator_m0s != NULL) {
        delete interpolator_m0s;
    }
}

void InversionBase::set_mesh(const Mesh& mesh0) {
    this->mesh = mesh0;
    this->Nm = this->mesh.n_elems();
    mesh.set_n_parameter(5);
    m0 = VectorXd::Constant(Nm, 1, 0);
    m.resize(Nm);
    m_min = VectorXd::Constant(Nm, 1, -1e6);
    m_max = VectorXd::Constant(Nm, 1, 1e6);
    m_ini = VectorXd::Constant(Nm, 1, 0);
    m0_s = VectorXd::Constant(Nm, 1, 0);

    this->init_matrices();
}

void InversionBase::set_S() {
    S_s.resize(Nm, Nm);  // Resizes the matrix to a rows x cols matrix and
    // initializes it to zero.
    S_s.reserve(Nm);
    // S_s.setZero();

    S_theta.resize(Nm, Nm);
    S_theta.reserve(Nm);
    // S_theta.setZero();

    S_phi.resize(Nm, Nm);
    S_phi.reserve(Nm);
    // S_phi.setZero();

    S_r.resize(Nm, Nm);
    S_r.reserve(Nm);
    // S_r.setZero();

    S_crg.resize(Nm, Nm);
    S_crg.reserve(Nm);

    // S_crg.setZero();

    double dr, dtheta, dphi;
    double r0, theta0, phi0;

    // mesh->leaf_cells[0]->get_size(dr,dtheta,dphi);
    // double min_side=dr;

    for (int i = 0; i < Nm; i++) {
        Cell* c = mesh.leaf_cells[i];
        c->get_size(dr, dtheta, dphi);
        c->get_center(r0, theta0, phi0);

        S_s.coeffRef(i, i) += 1;
        S_theta.coeffRef(i, i) += 1;
        S_phi.coeffRef(i, i) += 1;
        S_r.coeffRef(i, i) += 1;
        S_crg.coeffRef(i, i) += 1;
        // S_s.coeffRef(i, i) += 1;
        // S_theta.coeffRef(i, i) += r0 * dtheta;
        // S_phi.coeffRef(i, i) += r0 * sin(theta0) * dphi;
        // S_r.coeffRef(i, i) += dr;

        // min_side=(min_side>dr)?dr:min_side;
        // min_side=(min_side>(r0 * dtheta))?(r0 * dtheta):min_side;
        // min_side=(min_side>(r0 * sin(theta0) * dphi))?(r0 * sin(theta0) *
        // dphi):min_side;
    }
    // S_theta=S_theta/min_side;
    // S_phi=S_phi/min_side;
    // S_r=S_r/min_side;

    S_s.makeCompressed();
    S_theta.makeCompressed();
    S_phi.makeCompressed();
    S_r.makeCompressed();
    S_crg.makeCompressed();
    S_s.data().squeeze();
    S_theta.data().squeeze();
    S_phi.data().squeeze();
    S_r.data().squeeze();
    S_crg.data().squeeze();
}

void InversionBase::init_matrices() {
    // This is necesseary because cell number is changed after the mesh is
    // refines
    this->Nm = this->mesh.n_elems();
    // cout << Nm << endl;

    V.resize(Nm, Nm);
    V.reserve(Nm);
    this->set_Vmatrix();
    V.makeCompressed();
    V.data().squeeze();

    R_r.resize(Nm, Nm);
    R_r.reserve(Nm);
    // R_r.setZero();
    this->set_depth_weighting(this->depth_weighting_factor);
    R_r.makeCompressed();
    R_r.data().squeeze();

    this->set_S();

    D_s.resize(Nm, Nm);  // resize and initialize it to zeros
    D_s.reserve(Nm);
    // D_s.setZero();

    D_theta1.resize(Nm, Nm);
    D_theta1.reserve(Nm * 2);
    // D_theta1.setZero();

    D_phi1.resize(Nm, Nm);
    D_phi1.reserve(Nm * 2);
    // D_phi1.setZero();

    D_r1.resize(Nm, Nm);
    D_r1.reserve(Nm * 2);
    // D_r1.setZero();
    this->set_difference_matrix();
    D_s.makeCompressed();
    D_s.data().squeeze();
    D_theta1.makeCompressed();
    D_theta1.data().squeeze();
    D_phi1.makeCompressed();
    D_phi1.data().squeeze();
    D_r1.makeCompressed();
    D_r1.data().squeeze();
}

void InversionBase::set_difference_matrix() {
    D_s.setZero();
    D_theta1.setZero();
    D_phi1.setZero();
    D_r1.setZero();

    //#pragma omp parallel for
    for (int i = 0; i < Nm; i++) {
        Cell* c = mesh.leaf_cells[i];
        int id0 = c->get_id();
        double dr0, dtheta0, dphi0, dr1, dtheta1, dphi1;
        double r0c, theta0c, phi0c;
        double r1c, theta1c, phi1c;
        c->get_center(r0c, theta0c, phi0c);
        c->get_size(dr0, dtheta0, dphi0);

        D_s.coeffRef(id0, id0) += 1.0;

        Face* f = c->external_faces_theta[1];
        if (f->isleaf) {
            Cell* neigh = f->neigh_cells[1];
            if (f->neigh_cells[0] != NULL && f->neigh_cells[1] != NULL) {
                // if()
                neigh->get_size(dr1, dtheta1, dphi1);
                int id1 = neigh->get_id();
                assert(id0 != id1);
                D_theta1.coeffRef(id0, id0) +=
                    -1.0 / (0.5 * r0c * (dtheta0 + dtheta1));
                D_theta1.coeffRef(id0, id1) +=
                    1.0 / (0.5 * r0c * (dtheta0 + dtheta1));
            } else {
                f = c->external_faces_theta[0];
                if (f->isleaf) {
                    neigh = f->neigh_cells[0];
                    // if (neigh == NULL)
                    // {
                    //     neigh = f->neigh_cells[1];
                    // }
                    neigh->get_size(dr1, dtheta1, dphi1);
                    int id1 = neigh->get_id();
                    assert(id0 != id1);
                    D_theta1.coeffRef(id0, id0) +=
                        1.0 / (0.5 * r0c * (dtheta0 + dtheta1));
                    D_theta1.coeffRef(id0, id1) +=
                        -1.0 / (0.5 * r0c * (dtheta0 + dtheta1));
                } else {
                    Cell* neigh = f->child_faces[0]->neigh_cells[0];
                    neigh->get_size(dr1, dtheta1, dphi1);
                    D_theta1.coeffRef(id0, id0) +=
                        1.0 / (0.5 * r0c * (dtheta0 + dtheta1));
                    for (int k = 0; k < 4; k++) {
                        Cell* neigh = f->child_faces[k]->neigh_cells[0];
                        int id1 = neigh->get_id();
                        assert(id0 != id1);
                        D_theta1.coeffRef(id0, id1) +=
                            -1.0 / (0.5 * r0c * (dtheta0 + dtheta1)) * 0.25;
                    }
                }
            }
        } else {
            Cell* neigh = f->child_faces[0]->neigh_cells[1];
            neigh->get_size(dr1, dtheta1, dphi1);
            D_theta1.coeffRef(id0, id0) +=
                -1.0 / (0.5 * r0c * (dtheta0 + dtheta1));
            for (int k = 0; k < 4; k++) {
                Cell* neigh = f->child_faces[k]->neigh_cells[1];
                int id1 = neigh->get_id();
                assert(id0 != id1);
                D_theta1.coeffRef(id0, id1) +=
                    1.0 / (0.5 * r0c * (dtheta0 + dtheta1)) * 0.25;
            }
        }
        // y

        f = c->external_faces_phi[1];
        if (f->isleaf) {
            Cell* neigh = f->neigh_cells[1];
            if (f->neigh_cells[0] != NULL && f->neigh_cells[1] != NULL) {
                int id1 = neigh->get_id();
                neigh->get_size(dr1, dtheta1, dphi1);
                assert(id0 != id1);
                D_phi1.coeffRef(id0, id0) +=
                    -1.0 / (0.5 * r0c * sin(theta0c) * (dphi0 + dphi1));
                D_phi1.coeffRef(id0, id1) +=
                    1.0 / (0.5 * r0c * sin(theta0c) * (dphi0 + dphi1));
            } else {
                f = c->external_faces_phi[0];
                if (f->isleaf) {
                    neigh = f->neigh_cells[0];
                    // if (neigh == NULL)
                    // {
                    //     neigh = c->external_faces_phi[0]->neigh_cells[1];
                    // }
                    neigh->get_size(dr1, dtheta1, dphi1);
                    int id1 = neigh->get_id();
                    assert(id0 != id1);
                    D_phi1.coeffRef(id0, id0) +=
                        1.0 / (0.5 * r0c * sin(theta0c) * (dphi0 + dphi1));
                    D_phi1.coeffRef(id0, id1) +=
                        -1.0 / (0.5 * r0c * sin(theta0c) * (dphi0 + dphi1));
                } else {
                    Cell* neigh = f->child_faces[0]->neigh_cells[0];
                    neigh->get_size(dr1, dtheta1, dphi1);
                    D_phi1.coeffRef(id0, id0) +=
                        1.0 / (0.5 * r0c * sin(theta0c) * (dphi0 + dphi1));
                    for (int k = 0; k < 4; k++) {
                        Cell* neigh = f->child_faces[k]->neigh_cells[0];
                        int id1 = neigh->get_id();
                        assert(id0 != id1);
                        D_phi1.coeffRef(id0, id1) +=
                            -1.0 /
                            (0.5 * r0c * sin(theta0c) * (dphi0 + dphi1)) * 0.25;
                    }
                }
            }
        } else {
            Cell* neigh = f->child_faces[0]->neigh_cells[1];
            neigh->get_size(dr1, dtheta1, dphi1);
            D_phi1.coeffRef(id0, id0) +=
                -1.0 / (0.5 * r0c * sin(theta0c) * (dphi0 + dphi1));
            for (int k = 0; k < 4; k++) {
                Cell* neigh = f->child_faces[k]->neigh_cells[1];
                int id1 = neigh->get_id();
                assert(id0 != id1);
                D_phi1.coeffRef(id0, id1) +=
                    1.0 / (0.5 * r0c * sin(theta0c) * (dphi0 + dphi1)) * 0.25;
            }
        }

        // z
        f = c->external_faces_r[1];
        if (f->isleaf) {
            Cell* neigh = f->neigh_cells[1];
            if (f->neigh_cells[0] != NULL && f->neigh_cells[1] != NULL) {
                int id1 = neigh->get_id();
                neigh->get_size(dr1, dtheta1, dphi1);
                if (id0 == id1) {
                    cout << "c" << endl;
                    cout << "North-south direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        c->external_faces_theta[j]->display();
                    }
                    cout << "West-east direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        c->external_faces_phi[j]->display();
                    }
                    cout << "Radius direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        c->external_faces_r[j]->display();
                    }
                    cout << "-----------------------------" << endl;
                    cout << "neigh 0" << endl;
                    cout << "North-south direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        f->neigh_cells[0]->external_faces_theta[j]->display();
                    }
                    cout << "West-east direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        f->neigh_cells[0]->external_faces_phi[j]->display();
                    }
                    cout << "Radius direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        f->neigh_cells[0]->external_faces_r[j]->display();
                    }
                    cout << "-----------------------------" << endl;
                    cout << "neigh 1" << endl;
                    cout << "North-south direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        f->neigh_cells[1]->external_faces_theta[j]->display();
                    }
                    cout << "West-east direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        f->neigh_cells[1]->external_faces_phi[j]->display();
                    }
                    cout << "Radius direction" << endl;
                    for (int j = 0; j < 2; j++) {
                        f->neigh_cells[1]->external_faces_r[j]->display();
                    }
                }
                assert(id0 != id1);

                D_r1.coeffRef(id0, id0) += -1.0 / (0.5 * (dr0 + dr1));
                D_r1.coeffRef(id0, id1) += 1.0 / (0.5 * (dr0 + dr1));
            } else if (c->external_faces_r[0]->neigh_cells[0] != NULL) {
                f = c->external_faces_r[0];
                if (f->isleaf) {
                    neigh = f->neigh_cells[0];
                    neigh->get_size(dr1, dtheta1, dphi1);
                    int id1 = neigh->get_id();
                    assert(id0 != id1);
                    D_r1.coeffRef(id0, id0) += 1.0 / (0.5 * (dr0 + dr1));
                    D_r1.coeffRef(id0, id1) += -1.0 / (0.5 * (dr0 + dr1));
                } else {
                    Cell* neigh = f->child_faces[0]->neigh_cells[0];
                    neigh->get_size(dr1, dtheta1, dphi1);
                    D_r1.coeffRef(id0, id0) += 1.0 / (0.5 * (dr0 + dr1));
                    for (int k = 0; k < 4; k++) {
                        Cell* neigh = f->child_faces[k]->neigh_cells[0];
                        neigh->get_size(dr1, dtheta1, dphi1);
                        int id1 = neigh->get_id();
                        assert(id0 != id1);
                        D_r1.coeffRef(id0, id1) +=
                            -1.0 / (0.5 * (dr0 + dr1)) * 0.25;
                    }
                }
            }
        } else {
            Cell* neigh = f->child_faces[0]->neigh_cells[1];
            neigh->get_size(dr1, dtheta1, dphi1);
            D_r1.coeffRef(id0, id0) += -1.0 / (0.5 * (dr0 + dr1));
            for (int k = 0; k < 4; k++) {
                Cell* neigh = f->child_faces[k]->neigh_cells[1];
                neigh->get_size(dr1, dtheta1, dphi1);
                int id1 = neigh->get_id();
                assert(id0 != id1);
                D_r1.coeffRef(id0, id1) += 1.0 / (0.5 * (dr0 + dr1)) * 0.25;
            }
        }
    }
}

void InversionBase::show_differece_matrix(unsigned int direction) {
    if (direction == NORTH_SOUTH) {
        cout << D_theta1 << endl;
    } else if (direction == WEST_EAST) {
        cout << D_phi1 << endl;
    } else if (direction == RADIUS) {
        cout << D_r1 << endl;
    }
}

void InversionBase::set_Vmatrix() {
    V.setZero();
    // V.setIdentity();
    assert(N_obs != 0);
    assert(Nm != 0);
    double max_v = 0;
    double total_v = 0;
    for (int i = 0; i < mesh.n_elems(); i++) {
        const Tesseroid& t = mesh.get_elem(i);
        double v = t.get_volumn();
        V.coeffRef(i, i) += sqrt(v);
        // V.coeffRef(i, i) += 1;
    }
}

void InversionBase::set_Tmatrix() {
    T_r.resize(Nm, Nm);
    T_r.data().squeeze();
    T_theta.resize(Nm, Nm);
    T_theta.data().squeeze();
    T_phi.resize(Nm, Nm);
    T_phi.data().squeeze();
    VectorXd Drs_v = D_r1 * m0_s;
    VectorXd Dthetas_v = D_theta1 * m0_s;
    VectorXd Dphis_v = D_phi1 * m0_s;

    SMatrix Drs(Nm, Nm);
    SMatrix Dthetas(Nm, Nm);
    SMatrix Dphis(Nm, Nm);
    Drs.setZero();
    Dthetas.setZero();
    Dphis.setZero();
    for (int i = 0; i < Nm; i++) {
        Drs.coeffRef(i, i) = Drs_v(i);
        Dthetas.coeffRef(i, i) = Dthetas_v(i);
        Dphis.coeffRef(i, i) = Dphis_v(i);
    }
    T_r = V * (Dphis * D_theta1 - Dthetas * D_phi1);
    T_theta = V * (Drs * D_phi1 - Dphis * D_r1);
    T_phi = V * (Dthetas * D_r1 - Drs * D_theta1);

    T_r.makeCompressed();
    T_theta.makeCompressed();
    T_phi.makeCompressed();
}

void InversionBase::set_depth_weighting(double beta, int flag) {
    this->depth_weighting_factor = beta;
    R_r.setZero();
    assert(N_obs != 0);  // Nd>=N_obs>=0
    assert(Nm != 0);
    if (flag == 0) {
        double rp = 0;
        for (int i = 0; i < N_obs; i++) {
            Point p = ob(i);
            rp += p(0);
        }
        rp = rp / N_obs;
        //#pragma omp parallel for
        for (int i = 0; i < Nm; i++) {
            const Tesseroid& t = mesh.get_elem(i);
            double r = 0.5 * (t._r[0] + t._r[1]);
            double weight = 1.0 / pow(std::abs(rp - r), beta * 0.5);
            R_r.coeffRef(i, i) += weight;
        }
    }
    if (flag == 1) {
        assert(use_wavelet = false);
        for (int i = 0; i < Nm; i++) {
            double t = (G.col(i)).norm() / Nd;
            R_r.coeffRef(i, i) += sqrt(t);
        }
    }
}

void InversionBase::set_dobs(const VectorXd& d) {
    int n = field_flag.count();
    assert(Nd == d.size());

    this->dobs = d;
    // cout<<dobs.rows()<<endl;
    Wd.resize(Nd, Nd);
    Wd.setZero();
    for (int i = 0; i < Nd; i++) {
        Wd.coeffRef(i, i) += 1.0;
    }
}

void InversionBase::set_dobs(const VectorXd& d, double relative_error,
                             double a) {
    assert(relative_error < 0.7);
    int n = field_flag.count();
    assert(Nd == d.size());

    this->dobs = d;
    Wd.resize(Nd, Nd);
    Wd.setZero();
    for (int i = 0; i < Nd; i++) {
        Wd.coeffRef(i, i) += 1.0 / (relative_error * std::fabs(dobs(i)) + a);
    }
}
void InversionBase::set_dobs(const VectorXd& d, vector<double> relative_error,
                             vector<double> a) {
    int n_fields = field_flag.count();
    assert(relative_error.size() == n_fields && a.size() == n_fields);
    assert(Nd == d.size());
    this->dobs = d;
    Wd.resize(Nd, Nd);
    Wd.setZero();
    for (int i = 0; i < n_fields; i++) {
        for (int j = 0; j < N_obs; j++) {
            Wd.coeffRef(j + i * N_obs, j + i * N_obs) +=
                1.0 /
                (relative_error[i] * std::fabs(dobs(j + i * N_obs)) + a[i]);
        }
    }
}

void InversionBase::set_dobs2(const VectorXd& d, double a) {
    assert(a > 0);
    int n = field_flag.count();
    assert(Nd == d.size());

    this->dobs = d;
    Wd.resize(Nd, Nd);
    Wd.setZero();
    for (int i = 0; i < Nd; i++) {
        Wd.coeffRef(i, i) += 1.0 / a;
    }
}

void InversionBase::set_Wd(const VectorXd& sigma) {
    assert(sigma.rows() == Wd.cols());
    Wd.setZero();
    for (int i = 0; i < Nd; i++) {
        Wd.coeffRef(i, i) += 1.0 / sigma(i);
    }
}

void InversionBase::set_observation(const Observation& ob0) {
    this->ob = ob0;
    N_obs = ob.get_n_obs();
    int n_fields = field_flag.count();
    Nd = N_obs * n_fields;

    Wd.resize(Nd, Nd);
    Wd.reserve(Nd);
    Wd.setIdentity();
}

void InversionBase::set_weights_of_objectives(double as, double ar,
                                              double a_theta, double a_phi,
                                              double a_crg) {
    this->a_s = as;
    this->a_r = ar;
    this->a_theta = a_theta;
    this->a_phi = a_phi;
    this->a_crg = a_crg;
}

void InversionBase::set_reference_model(VectorXd& m_ref) { this->m0 = m_ref; }

void InversionBase::set_geometry_reference_model(VectorXd& m_ref) {
    this->m0_s = m_ref;
}
void InversionBase::set_m(VectorXd& m_) { this->m = m_; }

void InversionBase::set_m_ini(VectorXd& m_ini_) { this->m_ini = m_ini_; }

void InversionBase::set_min_max(VectorXd& m_min_, VectorXd& m_max_) {
    this->m_min = m_min_;
    this->m_max = m_max_;
}

void InversionBase::set_density_to_mesh() {
    assert(m.size() == Nm);
#pragma omp parallel for
    for (int i = 0; i < Nm; i++) {
        Cell* t = mesh.leaf_cells[i];
        t->set_parameter(m(i), 0);
        // t.set_density(log10(R_r.coeffRef(i,i)));
    }
}

void InversionBase::set_reference_model_to_mesh() {
    assert(m0.size() == Nm);
#pragma omp parallel for
    for (int i = 0; i < Nm; i++) {
        Cell* t = mesh.leaf_cells[i];
        assert(t->parameters.size() > 1);
        t->set_parameter(m0(i), 1);
        t->set_parameter(m0_s(i), 2);
        // t.set_density(log10(Z.coeffRef(i,i)));
    }
}

void InversionBase::set_min_max_to_mesh() {
    assert(m.size() == Nm);
    for (int i = 0; i < Nm; i++) {
        Cell* t = mesh.leaf_cells[i];
        t->set_parameter(m_min(i), 3);
        t->set_parameter(m_max(i), 4);
        // t.set_density(log10(Z.coeffRef(i,i)));
    }
}

void InversionBase::set_petrophysics_constraint(
    VectorXd& m_ref_other_para, function<double(double)> relation) {
    int ns = m_ref_other_para.size();
    assert(ns == Nm);
    for (int i = 0; i < Nm; i++) {
        // if (fabs(m_ref_other_para(i)) < 1e-9)
        // {
        //     this->m0(i) = 0;
        // }
        // else
        // {
        this->m0(i) = relation(m_ref_other_para(i));
        // }
    }
}

VectorXd InversionBase::get_predicted_field() {
    VectorXd d_pre;
    if (use_wavelet) {
        this->G_vec_mul(this->m, d_pre);
    } else {
        d_pre = (this->G) * (this->m);
    }
    return d_pre;
}

void InversionBase::update_S_crg() {
    S_crg.setZero();
    for (int i_CELL = 0; i_CELL < Nm; i_CELL++) {
        double rc, thetac, phic;
        double depthc;
        mesh.leaf_cells[i_CELL]->get_center(rc, thetac, phic);
        double latc = 90.0 - thetac * 180.0 / GS::PI;
        double lonc = phic * 180.0 / GS::PI - 180.0;

        depthc = (reference_surface - rc) / 1000.0;
        if ((depthc < constraint_depth[1] ||
             std::abs(depthc - constraint_depth[1]) < 1e-10) &&
            (depthc > constraint_depth[0] ||
             std::abs(depthc - constraint_depth[0]) < 1e-10) &&
            (latc < constraint_lat[1] ||
             std::abs(latc - constraint_lat[1]) < 1e-10) &&
            (latc > constraint_lat[0] ||
             std::abs(latc - constraint_lat[0]) < 1e-10) &&
            (lonc < constraint_lon[1] ||
             std::abs(lonc - constraint_lon[1]) < 1e-10) &&
            (lonc > constraint_lon[0] ||
             std::abs(lonc - constraint_lon[0]) < 1e-10)) {
            S_crg.coeffRef(i_CELL, i_CELL) += 1;
        }
    }
}

void InversionBase::output_predicted_data(string out_name) {
    VectorXd d_pre = this->get_predicted_field();
    this->out_data(d_pre, out_name);
}

void InversionBase::output_obs_data(string out_name) {
    // cout<<this->dobs.rows()<<endl;
    this->out_data(this->dobs, out_name);
}

void InversionBase::out_data(const VectorXd& d, string out_name) {
    vector<string> strs = {"V",          "g_r",      "g_theta", "g_phi",
                           "T_rr",       "T_rtheta", "T_rphi",  "T_thetatheta",
                           "T_thetaphi", "T_phiphi"};
    vector<unsigned int> field_label;
    for (unsigned int i = 0; i < strs.size(); i++) {
        if (field_flag[i]) {
            field_label.push_back(i);
        }
    }

    int n_com = field_flag.count();  // number of used components
    int n_ob = ob.get_n_obs();
    int nd = d.rows();
    // cout<<"xx  "<<d.rows()<<endl;

    // cout<<nd<<endl<<n_ob<<endl<<n_com<<endl;

    assert(nd % n_ob == 0);
    assert(nd / n_ob == n_com);

    for (int j = 0; j < n_com; j++) {
        string file_name = out_name + "_" + strs[field_label[j]];
        ofstream out_s(file_name);
        for (int i = 0; i < n_ob; i++) {
            // get_n_obs() obtains number of observation points, data size =
            // n_com*
            const Point& p = ob(i);
            out_s << fixed;
            out_s << setw(25) << setprecision(7) << left
                  << 90.0 - p(1) * 180.0 / GS::PI << setw(25) << left
                  << p(2) * 180.0 / GS::PI - 180.0;
            out_s << scientific;
            out_s << setw(30) << setprecision(15) << left << d(i + j * n_ob)
                  << endl;
        }
    }
}

void InversionBase::cumulative_sensitivity() {
    VectorXd sensitivity(Nm);
    cout << "Calculating cumulative sensitivity ..." << endl;
    for (int i = 0; i < Nm; i++) {
        Cell* t = mesh.leaf_cells[i];
        sensitivity(i) = 0;
        for (int j = 0; j < Nd; j++) {
            sensitivity(i) += std::abs(G(j, i));
        }
        double size = t->get_volumn();
        sensitivity(i) = sensitivity(i) / size;
    }
    double max_s = sensitivity.maxCoeff();
    double min_s = sensitivity.minCoeff();
    sensitivity =
        (sensitivity - VectorXd::Constant(Nm, min_s)) / (max_s - min_s);

    // occupy the first parameter of the cells temporarily
    for (int i = 0; i < Nm; i++) {
        Cell* t = mesh.leaf_cells[i];
        t->set_parameter(sensitivity(i), 0);
    }
    vector<string> parameter_name = {"sensitivity"};
    this->mesh.out_model_netcdf("cumulative_sensitivity.nc", 0, "sensitivity",
                                "");
    this->mesh.out_model_vtk(string("cumulative_sensitivity.vtk"), 1,
                             parameter_name);
    cout << "The values of cumulative sensitivity in each cell have been "
            "written into file cumulative_sensitivity.nc and "
            "cumulative_sensitivity.vtk"
         << endl;
    // store density values into the first parameter again
    this->set_density_to_mesh();
}

void InversionBase::result2txt(string filename) {
    this->set_density_to_mesh();
    mesh.out_model_txt(filename + string(".txt"));
}

void InversionBase::result2vtk(string filename) {
    this->set_density_to_mesh();
    this->set_reference_model_to_mesh();
    vector<string> parameter_name = {"model", "m0", "m0s"};
    mesh.out_model_vtk(filename + string(".vtk"), 3, parameter_name);
}
void InversionBase::result2netcdf(string filename) {
    this->set_density_to_mesh();
    this->mesh.out_model_netcdf(filename + string(".nc"));
}

void InversionBase::create_crg_model_from_data(string filename, int lat_size,
                                               int lon_size, int dep_size,
                                               string data_order,
                                               int fast_dimension,
                                               double reference_level) {
    this->use_cross_gradient_constraint = true;
    if (this->interpolator_m0s != NULL) {
        delete interpolator_m0s;
        interpolator_m0s = NULL;
    }

    ifstream input_file;
    input_file.open(filename);
    assert(input_file.good());
    cout << "Read cross-gradient constraint model from " << filename << endl;

    vector<double> grid_x;  // latitude, unit: degree
    vector<double> grid_y;  // longitude, unit: degree
    vector<double> grid_z;  // depth, uniy: kilometer

    grid_x.resize(lat_size);
    grid_y.resize(lon_size);
    grid_z.resize(dep_size);

    // the size of the grid in each dimension
    array<int, 3> grid_sizes;
    grid_sizes[0] = grid_x.size();
    grid_sizes[1] = grid_y.size();
    grid_sizes[2] = grid_z.size();

    int num_elements = grid_sizes[0] * grid_sizes[1] * grid_sizes[2];
    std::vector<double> f_values(num_elements);

    double lat, lon, dep, val;  // each column

    if (fast_dimension == 0) {
        for (int k = 0; k < grid_z.size(); k++) {
            for (int j = 0; j < grid_y.size(); j++) {
                for (int i = 0; i < grid_x.size(); i++) {
                    if (data_order == "yxz") {
                        input_file >> lon >> lat >> dep >> val;
                    } else if (data_order == "xyz") {
                        input_file >> lat >> lon >> dep >> val;
                    } else if (data_order == "zxy") {
                        input_file >> dep >> lat >> lon >> val;
                    } else if (data_order == "zyx") {
                        input_file >> dep >> lon >> lat >> val;
                    } else if (data_order == "yzx") {
                        input_file >> lon >> dep >> lat >> val;
                    } else if (data_order == "xzy") {
                        input_file >> lat >> dep >> lon >> val;
                    } else {
                        cout << "data_order should be one of the following:"
                             << endl;
                        cout << "xyz yxz zxy zyx xzy yzx" << endl;
                        std::abort();
                    }

                    if (k == 0 && j == 0) {
                        grid_x[i] = lat;  // latitude
                    }

                    f_values[i * grid_sizes[1] * grid_sizes[2] +
                             j * grid_sizes[2] + k] = val;
                }
                if (k == 0) {
                    grid_y[j] = lon;  // longitude
                }
            }
            grid_z[k] = dep;  // depth
        }
    } else if (fast_dimension == 1) {
        for (int k = 0; k < grid_z.size(); k++) {
            for (int i = 0; i < grid_x.size(); i++) {
                for (int j = 0; j < grid_y.size(); j++) {
                    if (data_order == "yxz") {
                        input_file >> lon >> lat >> dep >> val;
                    } else if (data_order == "xyz") {
                        input_file >> lat >> lon >> dep >> val;
                    } else if (data_order == "zxy") {
                        input_file >> dep >> lat >> lon >> val;
                    } else if (data_order == "zyx") {
                        input_file >> dep >> lon >> lat >> val;
                    } else if (data_order == "yzx") {
                        input_file >> lon >> dep >> lat >> val;
                    } else if (data_order == "xzy") {
                        input_file >> lat >> dep >> lon >> val;
                    } else {
                        cout << "data_order should be one of the following:"
                             << endl;
                        cout << "xyz yxz zxy zyx xzy yzx" << endl;
                        std::abort();
                    }

                    if (k == 0 && i == 0) {
                        grid_y[j] = lon;  // longitude
                    }

                    f_values[i * grid_sizes[1] * grid_sizes[2] +
                             j * grid_sizes[2] + k] = val;
                }
                if (k == 0) {
                    grid_x[i] = lat;  // latitude
                }
            }
            grid_z[k] = dep;  // depth
        }
    } else {
        cout << "Parameter fast_dimension"
             << " should be 0 or 1" << endl;
    }
    input_file.close();
    input_file.clear();

    // construct the grid in each dimension.
    // note that we will pass in a sequence of iterators pointing to the
    // beginning of each grid
    std::vector<std::vector<double>::iterator> grid_iter_list;
    grid_iter_list.push_back(grid_x.begin());
    grid_iter_list.push_back(grid_y.begin());
    grid_iter_list.push_back(grid_z.begin());

    // construct the interpolator. the last two arguments are pointers to the
    // underlying data
    this->interpolator_m0s = new InterpMultilinear<3, double>(
        grid_iter_list.begin(), grid_sizes.begin(), f_values.data(),
        f_values.data() + num_elements);

    Mesh mesh_crg;
    double lat_space =
        (grid_x[grid_x.size() - 1] - grid_x[0]) / (grid_x.size() - 1);
    double lon_space =
        (grid_y[grid_y.size() - 1] - grid_y[0]) / (grid_y.size() - 1);
    double r_space =
        (grid_z[grid_z.size() - 1] - grid_z[0]) / (grid_z.size() - 1);  // km
    double lat_model[2] = {grid_x[0] - 0.5 * lat_space,
                           grid_x[grid_x.size() - 1] + 0.5 * lat_space};
    double lon_model[2] = {grid_y[0] - 0.5 * lon_space,
                           grid_y[grid_y.size() - 1] + 0.5 * lon_space};
    double depth_model[2] = {
        (grid_z[0] - 0.5 * r_space) * 1000,
        (grid_z[grid_z.size() - 1] + 0.5 * r_space) * 1000};
    mesh_crg.generate_regular_geographic_mesh(lat_model, lat_size, lon_model,
                                              lon_size, depth_model, dep_size,
                                              reference_level);

    for (int i = 0; i < mesh_crg.n_elems(); i++) {
        Cell* c = mesh_crg.leaf_cells[i];
        double rc, thetac, phic;
        c->get_center(rc, thetac, phic);

        double latc = 90.0 - thetac * 180.0 / GS::PI;
        double lonc = phic * 180.0 / GS::PI - 180.0;
        double depthc = (reference_surface - rc) / 1000.0;

        array<double, 3> args = {latc, lonc, depthc};
        double val = (*interpolator_m0s).interp(args.begin());
        // cout<<"#########"<<c->parameters.size()<<endl;
        c->set_parameter(val);
    }
    mesh_crg.out_model_netcdf(string("crg_model.nc"), 0, "crg", "");
    mesh_crg.out_model_vtk("crg_model.vtk", 1, vector<string>(1, "crg"));

    // interpolator_m0sc
}

void InversionBase::create_ref_model_from_data(string filename, int lat_size,
                                               int lon_size, int dep_size,
                                               string data_order,
                                               int fast_dimension,
                                               double reference_level) {
    this->use_petrophysical_constraint = true;
    if (this->interpolator_m0 != NULL) {
        delete interpolator_m0;
        interpolator_m0 = NULL;
    }

    ifstream input_file;
    input_file.open(filename);
    assert(input_file.good());
    cout << "Read converted density model from " << filename << endl;

    vector<double> grid_x;  // latitude
    vector<double> grid_y;  // longitude
    vector<double> grid_z;  // depth

    grid_x.resize(lat_size);
    grid_y.resize(lon_size);
    grid_z.resize(dep_size);

    // the size of the grid in each dimension
    array<int, 3> grid_sizes;
    grid_sizes[0] = grid_x.size();
    grid_sizes[1] = grid_y.size();
    grid_sizes[2] = grid_z.size();

    int num_elements = grid_sizes[0] * grid_sizes[1] * grid_sizes[2];
    std::vector<double> f_values(num_elements);

    double lat, lon, dep, val;  // each column

    if (fast_dimension == 0) {
        for (int k = 0; k < grid_z.size(); k++) {
            for (int j = 0; j < grid_y.size(); j++) {
                for (int i = 0; i < grid_x.size(); i++) {
                    if (data_order == "yxz") {
                        input_file >> lon >> lat >> dep >> val;
                    } else if (data_order == "xyz") {
                        input_file >> lat >> lon >> dep >> val;
                    } else if (data_order == "zxy") {
                        input_file >> dep >> lat >> lon >> val;
                    } else if (data_order == "zyx") {
                        input_file >> dep >> lon >> lat >> val;
                    } else if (data_order == "yzx") {
                        input_file >> lon >> dep >> lat >> val;
                    } else if (data_order == "xzy") {
                        input_file >> lat >> dep >> lon >> val;
                    } else {
                        cout << "data_order should be one of the following:"
                             << endl;
                        cout << "xyz yxz zxy zyx xzy yzx" << endl;
                        std::abort();
                    }

                    if (k == 0 && j == 0) {
                        grid_x[i] = lat;  // latitude
                    }

                    f_values[i * grid_sizes[1] * grid_sizes[2] +
                             j * grid_sizes[2] + k] = val;
                }
                if (k == 0) {
                    grid_y[j] = lon;  // longitude
                }
            }
            grid_z[k] = dep;  // depth
        }
    } else if (fast_dimension == 1) {
        for (int k = 0; k < grid_z.size(); k++) {
            for (int i = 0; i < grid_x.size(); i++) {
                for (int j = 0; j < grid_y.size(); j++) {
                    if (data_order == "yxz") {
                        input_file >> lon >> lat >> dep >> val;
                    } else if (data_order == "xyz") {
                        input_file >> lat >> lon >> dep >> val;
                    } else if (data_order == "zxy") {
                        input_file >> dep >> lat >> lon >> val;
                    } else if (data_order == "zyx") {
                        input_file >> dep >> lon >> lat >> val;
                    } else if (data_order == "yzx") {
                        input_file >> lon >> dep >> lat >> val;
                    } else if (data_order == "xzy") {
                        input_file >> lat >> dep >> lon >> val;
                    } else {
                        cout << "data_order should be one of the following:"
                             << endl;
                        cout << "xyz yxz zxy zyx xzy yzx" << endl;
                        std::abort();
                    }

                    if (k == 0 && i == 0) {
                        grid_y[j] = lon;  // longitude
                    }

                    f_values[i * grid_sizes[1] * grid_sizes[2] +
                             j * grid_sizes[2] + k] = val;
                }
                if (k == 0) {
                    grid_x[i] = lat;  // latitude
                }
            }
            grid_z[k] = dep;  // depth
        }
    } else {
        cout << "Parameter fast_dimension"
             << " should be 0 or 1" << endl;
    }
    input_file.close();
    input_file.clear();

    // construct the grid in each dimension.
    // note that we will pass in a sequence of iterators pointing to the
    // beginning of each grid
    std::vector<std::vector<double>::iterator> grid_iter_list;
    grid_iter_list.push_back(grid_x.begin());
    grid_iter_list.push_back(grid_y.begin());
    grid_iter_list.push_back(grid_z.begin());

    // construct the interpolator. the last two arguments are pointers to the
    // underlying data
    this->interpolator_m0 = new InterpMultilinear<3, double>(
        grid_iter_list.begin(), grid_sizes.begin(), f_values.data(),
        f_values.data() + num_elements);

    Mesh mesh_ref;
    double lat_space =
        (grid_x[grid_x.size() - 1] - grid_x[0]) / (grid_x.size() - 1);
    double lon_space =
        (grid_y[grid_y.size() - 1] - grid_y[0]) / (grid_y.size() - 1);
    double r_space =
        (grid_z[grid_z.size() - 1] - grid_z[0]) / (grid_z.size() - 1);  // km
    double lat_model[2] = {grid_x[0] - 0.5 * lat_space,
                           grid_x[grid_x.size() - 1] + 0.5 * lat_space};
    double lon_model[2] = {grid_y[0] - 0.5 * lon_space,
                           grid_y[grid_y.size() - 1] + 0.5 * lon_space};
    double depth_model[2] = {
        (grid_z[0] - 0.5 * r_space) * 1000,
        (grid_z[grid_z.size() - 1] + 0.5 * r_space) * 1000};
    mesh_ref.generate_regular_geographic_mesh(lat_model, lat_size, lon_model,
                                              lon_size, depth_model, dep_size,
                                              reference_level);

    for (int i = 0; i < mesh_ref.n_elems(); i++) {
        Cell* c = mesh_ref.leaf_cells[i];
        double rc, thetac, phic;
        c->get_center(rc, thetac, phic);

        double latc = 90.0 - thetac * 180.0 / GS::PI;
        double lonc = phic * 180.0 / GS::PI - 180.0;
        double depthc = (reference_surface - rc) / 1000.0;

        array<double, 3> args = {latc, lonc, depthc};
        double val = (*interpolator_m0).interp(args.begin());
        c->set_parameter(val);
    }
    mesh_ref.out_model_netcdf(string("ref_model.nc"), 0, "ref", "");
    mesh_ref.out_model_vtk("ref_model.vtk", 1,
                           vector<string>(1, "reference_model"));
}
