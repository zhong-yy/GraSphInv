#include "GravityField.h"

GravityField::GravityField(int order_) {
    this->order = order_;
    this->unit_cube_gauss = new Gauss3D(this->order);
    this->unit_gauss_2D = new Gauss2D(this->order);
    this->unit_gauss_1D = new Gauss1D(this->order);
}

GravityField::~GravityField() {
    if (this->unit_cube_gauss != NULL) {
        delete this->unit_cube_gauss;
        this->unit_cube_gauss = NULL;
    }
    if (this->unit_gauss_2D != NULL) {
        delete this->unit_gauss_2D;
        this->unit_gauss_2D = NULL;
    }
    if (this->unit_gauss_1D != NULL) {
        delete this->unit_gauss_1D;
        this->unit_gauss_1D = NULL;
    }
}

void GravityField::g_dg_tg_own(const Point &x, Tesseroid *T, double &g,
                               Matrix31D &dg, Matrix3D &tg) {
    int size = 1;
    std::vector<Tesseroid *> Ts(size);
    Ts[0] = T;

    this->g_dg_tg_own_for_tesseroids(x, Ts, g, dg, tg);
}

void GravityField::g_dg_tg_own_for_tesseroids(const Point &x,
                                              std::vector<Tesseroid *> &Ts,
                                              double &g, Matrix31D &dg,
                                              Matrix3D &tg) {
    assert(Ts.size() > 0);
    int size = Ts.size();
    assert(size > 0);
    double r = x(0);
    double theta = x(1);
    double phi = x(2);

    std::vector<double> gs(size);
    std::vector<Matrix31D> dgs(size);
    std::vector<Matrix3D> tgs(size);

    g = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            tg(i, j) = 0.;  // VERY IMPORTANT, MUST BE SET TO 0.!
    for (int i = 0; i < 3; i++) dg(i) = 0.;

    const std::vector<Point> &q = unit_cube_gauss->get_points();
    const std::vector<double> &w = unit_cube_gauss->get_weights();
    assert(q.size() == w.size());
    assert(q.size() > 0);

    // loop each Tesseroid
    // #pragma omp parallel for
    for (unsigned int tess = 0; tess < size; tess++) {
        Tesseroid &T = (*Ts[tess]);
        double r1 = T._r[0];
        double r2 = T._r[1];
        double theta1 = T._theta[0];
        double theta2 = T._theta[1];
        double phi1 = T._phi[0];
        double phi2 = T._phi[1];
        double density = T._density;
        double volume = (r2 - r1) * (theta2 - theta1) * (phi2 - phi1);
        assert(volume > 0.);
        const double factor = density * G0;

        double temp_g = 0.;
        Matrix31D temp_dg;
        Matrix3D temp_tg;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                temp_tg(i, j) = 0.;  // VERY IMPORTANT, MUST BE SET TO 0.!
        for (int i = 0; i < 3; i++) temp_dg(i) = 0.;

        for (unsigned int p = 0; p < w.size(); p++) {
            const double WJ = w[p] * volume * 0.125;
            double rprime =
                0.5 * (1.0 - q[p](0)) * r1 + 0.5 * (1.0 + q[p](0)) * r2;
            double thetaprime =
                0.5 * (1.0 - q[p](1)) * theta1 + 0.5 * (1.0 + q[p](1)) * theta2;
            double phiprime =
                0.5 * (1.0 - q[p](2)) * phi1 + 0.5 * (1.0 + q[p](2)) * phi2;
            double K = rprime * rprime * std::sin(thetaprime);
            double cosPHI = std::sin(theta) * std::sin(thetaprime) *
                                std::cos(phi - phiprime) +
                            std::cos(theta) * std::cos(thetaprime);
            double delta1 = rprime * cosPHI - r;
            double delta2 = -1.0 * rprime *
                            (std::sin(theta) * std::cos(thetaprime) -
                             std::cos(theta) * std::sin(thetaprime) *
                                 std::cos(phi - phiprime));
            double delta3 =
                rprime * std::sin(thetaprime) * std::sin(phiprime - phi);
            double delta[3] = {delta1, delta2, delta3};
            double l =
                std::sqrt(delta1 * delta1 + delta2 * delta2 + delta3 * delta3);

            temp_g += WJ * K * (1.0 / l) * factor;
            for (unsigned int i = 0; i < 3; i++)
                temp_dg(i) += factor * WJ * K * delta[i] / (l * l * l);
            for (unsigned int i = 0; i < 3; i++) {
                for (unsigned int j = 0; j < 3; j++) {
                    double deltaij = 0.;
                    if (i == j) deltaij = 1.0;
                    double l3 = 1.0 / (l * l * l);
                    double l2 = 1.0 / (l * l);
                    double deltaideltaj = delta[i] * delta[j];
                    double value = factor * WJ * K * l3 *
                                   (3 * deltaideltaj * l2 - deltaij);
                    temp_tg(i, j) += value;
                }
            }
        }
        gs[tess] = temp_g;
        dgs[tess] = temp_dg;
        tgs[tess] = temp_tg;
    }
    for (unsigned int tess = 0; tess < size; tess++) {
        g = g + gs[tess];
        dg = dg + dgs[tess];
        tg = tg + tgs[tess];
    }

    return;
}

void GravityField::field_for_a_tesseroid(const Point &x, Tesseroid *tess,
                                         double rho, std::vector<double> &field,
                                         std::bitset<10> flag) {
    // assert(flag.any());
    std::vector<int> fields_needed;
    const unsigned count = flag.count();
    for (int i = 0; i < flag.size(); i++) {
        // bool status = flag[i];
        if (flag[i]) {
            fields_needed.push_back(i);
        }
    }
    // assert(fields_needed.size() == count);

    double (GravityField::*kernel[10])(const double &, const double *) = {
        &GravityField::kernel_V,          &GravityField::kernel_gr,
        &GravityField::kernel_gtheta,     &GravityField::kernel_gphi,
        &GravityField::kernel_T_rr,       &GravityField::kernel_T_rtheta,
        &GravityField::kernel_T_rphi,     &GravityField::kernel_T_thetatheta,
        &GravityField::kernel_T_thetaphi, &GravityField::kernel_T_phiphi};

    // flag: 0 potential, 1 gr, 2 g_theta, 3 g_phi, 4 T_rr, 5 T_rtheta, 6
    // T_rphi, 7 T_thetatheta, 8 T_thetaphi, 9 T_phiphi
    field.clear();
    field.resize(10);
    for (int i = 0; i < 10; i++) {
        field[i] = 0.0;
    }

    // computation point
    const double r = x(0);
    const double theta = x(1);
    const double phi = x(2);

    // source tesseroid
    Tesseroid &T = (*tess);
    const double r1 = T._r[0];
    const double r2 = T._r[1];
    const double theta1 = T._theta[0];
    const double theta2 = T._theta[1];
    const double phi1 = T._phi[0];
    const double phi2 = T._phi[1];
    const double density = rho;
    const double volume = (r2 - r1) * (theta2 - theta1) *
                          (phi2 - phi1);  // not real volume in geometric sense
    assert(volume > 0.);

    const std::vector<Point> &q = unit_cube_gauss->get_points();
    const std::vector<double> &w = unit_cube_gauss->get_weights();
    assert(q.size() == w.size());
    assert(q.size() > 0);

    const double factor = density * G0 * volume * 0.125;
    for (unsigned int p = 0; p < w.size(); p++) {
        const double WJ = w[p];
        double rprime = 0.5 * (1.0 - q[p](0)) * r1 + 0.5 * (1.0 + q[p](0)) * r2;
        double thetaprime =
            0.5 * (1.0 - q[p](1)) * theta1 + 0.5 * (1.0 + q[p](1)) * theta2;
        double phiprime =
            0.5 * (1.0 - q[p](2)) * phi1 + 0.5 * (1.0 + q[p](2)) * phi2;
        double K = rprime * rprime * std::sin(thetaprime);
        double cosPHI =
            std::sin(theta) * std::sin(thetaprime) * std::cos(phi - phiprime) +
            std::cos(theta) * std::cos(thetaprime);
        double delta1 = rprime * cosPHI - r;
        double delta2 =
            -1.0 * rprime *
            (std::sin(theta) * std::cos(thetaprime) -
             std::cos(theta) * std::sin(thetaprime) * std::cos(phi - phiprime));
        double delta3 =
            rprime * std::sin(thetaprime) * std::sin(phiprime - phi);
        double delta[3] = {delta1, delta2, delta3};
        double l =
            std::sqrt(delta1 * delta1 + delta2 * delta2 + delta3 * delta3);
        for (int i = 0; i < count; i++) {
            int index = fields_needed[i];
            field[index] += factor * WJ * K * (this->*kernel[index])(l, delta);
        }
    }

    for (int i = 1; i < 4; i++) {
        field[i] = field[i] * GS::SI2mGal;  //
    }

    for (int i = 4; i < 10; i++) {
        field[i] = field[i] * GS::SI2Eotvos;
    }
}

void GravityField::field_for_a_tesseroid_surf(const Point &x, Tesseroid *tess,
                                              double rho,
                                              std::vector<double> &field,
                                              std::bitset<10> flag) {
    // assert(flag.any());
    std::vector<int> fields_needed;
    const unsigned count = flag.count();
    for (int i = 0; i < flag.size(); i++) {
        // bool status = flag[i];
        if (flag[i]) {
            fields_needed.push_back(i);
        }
    }
    // assert(fields_needed.size() == count);

    // double (GravityField::*kernel[10])(const double &, const double *)
    // =
    //     {&GravityField::kernel_V,
    //      &GravityField::kernel_gr,
    //      &GravityField::kernel_gtheta,
    //      &GravityField::kernel_gphi,
    //      &GravityField::kernel_T_rr,
    //      &GravityField::kernel_T_rtheta,
    //      &GravityField::kernel_T_rphi,
    //      &GravityField::kernel_T_thetatheta,
    //      &GravityField::kernel_T_thetaphi,
    //      &GravityField::kernel_T_phiphi};
    double (GravityField::*kernelS[10])(const double &, const double &,
                                        const double &, const double &,
                                        const double &, const double &, int) = {
        &GravityField::kernel_V_S,         &GravityField::kernel_gr_S,
        &GravityField::kernel_gtheta_S,    &GravityField::kernel_gphi_S,
        &GravityField::kernel_Trr_S,       &GravityField::kernel_Trtheta_S,
        &GravityField::kernel_Trphi_S,     &GravityField::kernel_Tthetatheta_S,
        &GravityField::kernel_Tthetaphi_S, &GravityField::kernel_Tphiphi_S};
    // flag: 0 potential, 1 gr, 2 g_theta, 3 g_phi, 4 T_rr, 5 T_rtheta, 6
    // T_rphi, 7 T_thetatheta, 8 T_thetaphi, 9 T_phiphi
    field.clear();
    field.resize(10);
    for (int i = 0; i < 10; i++) {
        field[i] = 0.0;
    }

    // computation point
    const double ro = x(0);
    const double thetao = x(1);
    const double phio = x(2);

    // source tesseroid
    Tesseroid &T = (*tess);
    const double r1 = T._r[0];
    const double r2 = T._r[1];
    const double theta1 = T._theta[0];
    const double theta2 = T._theta[1];
    const double phi1 = T._phi[0];
    const double phi2 = T._phi[1];
    const double density = rho;
    // const double volume = (r2 - r1) * (theta2 - theta1) * (phi2 - phi1);
    // //not real volume in geometric sense assert(volume > 0.);

    const std::vector<double> &q1 = unit_gauss_1D->get_points();
    const std::vector<double> &w1 = unit_gauss_1D->get_weights();
    const std::vector<Point> &q2 = unit_gauss_2D->get_points();
    const std::vector<double> &w2 = unit_gauss_2D->get_weights();
    assert(q2.size() == w2.size());
    assert(q2.size() > 0);

    const double A = 0.25 * (theta2 - theta1) * (phi2 - phi1);
    const double B = 0.5 * (phi2 - phi1);
    const double C = 0.5 * (theta2 - theta1);
    double factor[10] = {0.5 * density * G0, -density * G0, -density * G0,
                         -density * G0,      density * G0,  density * G0,
                         density * G0,       density * G0,  density * G0,
                         density * G0};
    // double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    for (unsigned int p = 0; p < w2.size(); p++) {
        double theta =
            0.5 * (1.0 - q2[p](0)) * theta1 + 0.5 * (1.0 + q2[p](0)) * theta2;
        double phi =
            0.5 * (1.0 - q2[p](1)) * phi1 + 0.5 * (1.0 + q2[p](1)) * phi2;
        // if (flag[0])
        // {
        //   field[0] += factor[index] * A * w2[p] * (kernel_V_S(r2, theta, phi,
        //   ro, thetao, phio, 1) - kernel_V_S(r1, theta, phi, ro, thetao, phio,
        //   1));
        // }
        for (int i = 0; i < count; i++) {
            int index = fields_needed[i];
            field[index] +=
                factor[index] * A * w2[p] *
                ((this->*kernelS[index])(r2, theta, phi, ro, thetao, phio, 1) -
                 (this->*kernelS[index])(r1, theta, phi, ro, thetao, phio, 1));
        }
    }
    for (unsigned int p = 0; p < w1.size(); p++) {
        double phi = 0.5 * (1.0 - q1[p]) * phi1 + 0.5 * (1.0 + q1[p]) * phi2;
        double theta =
            0.5 * (1.0 - q1[p]) * theta1 + 0.5 * (1.0 + q1[p]) * theta2;

        for (int i = 0; i < count; i++) {
            int index = fields_needed[i];
            // field[index] += factor[index] * A * w2[p] *
            // ((this->*kernelS[index])(r2, theta, phi, ro, thetao, phio, 1) -
            // (this->*kernelS[index])(r1, theta, phi, ro, thetao, phio, 1));

            field[index] += factor[index] * B * w1[p] *
                            (((this->*kernelS[index])(r2, theta2, phi, ro,
                                                      thetao, phio, 2) -
                              (this->*kernelS[index])(r1, theta2, phi, ro,
                                                      thetao, phio, 2)) -
                             ((this->*kernelS[index])(r2, theta1, phi, ro,
                                                      thetao, phio, 2) -
                              (this->*kernelS[index])(r1, theta1, phi, ro,
                                                      thetao, phio, 2)));

            field[index] += factor[index] * C * w1[p] *
                            (((this->*kernelS[index])(r2, theta, phi2, ro,
                                                      thetao, phio, 3) -
                              (this->*kernelS[index])(r1, theta, phi2, ro,
                                                      thetao, phio, 3)) -
                             ((this->*kernelS[index])(r2, theta, phi1, ro,
                                                      thetao, phio, 3) -
                              (this->*kernelS[index])(r1, theta, phi1, ro,
                                                      thetao, phio, 3)));
        }
    }

    for (int i = 1; i < 4; i++) {
        field[i] = field[i] * GS::SI2mGal;  //
    }

    for (int i = 4; i < 10; i++) {
        field[i] = field[i] * GS::SI2Eotvos;
    }
}
