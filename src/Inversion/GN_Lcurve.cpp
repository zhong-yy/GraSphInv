#include "GN_Lcurve.h"
void GN_Lcurve::invert() {
  // set reference model using interpolation
  if (use_petrophysical_constraint == true) {
    for (int i = 0; i < Nm; i++) {
      Cell* c = this->mesh.leaf_cells[i];
      double rc, thetac, phic;
      c->get_center(rc, thetac, phic);

      double latc = 90.0 - thetac * 180.0 / GS::PI;
      double lonc = phic * 180.0 / GS::PI - 180.0;
      double depthc = (reference_surface - rc) / 1000.0;

      array<double, 3> args = {latc, lonc, depthc};
      double val = (*interpolator_m0).interp(args.begin());
      c->set_parameter(val, 1);
    }
    this->m0.resize(Nm);
    mesh.get_model_parameter_from_mesh(m0, 1);
  }

  if (use_cross_gradient_constraint == true) {
    cout << "Cross gradient constraint is used" << endl;
    for (int i = 0; i < Nm; i++) {
      Cell* c = this->mesh.leaf_cells[i];
      double rc, thetac, phic;
      c->get_center(rc, thetac, phic);

      double latc = 90.0 - thetac * 180.0 / GS::PI;
      double lonc = phic * 180.0 / GS::PI - 180.0;
      double depthc = (reference_surface - rc) / 1000.0;

      array<double, 3> args = {latc, lonc, depthc};
      double val = (*interpolator_m0s).interp(args.begin());
      c->set_parameter(val, 2);
    }
    this->m0_s.resize(Nm);
    mesh.get_model_parameter_from_mesh(m0_s, 2);
    this->set_Tmatrix();
    // update_S_crg();
  }

  SMatrix Ws = a_s * S_s * V * D_s * R_r;
  SMatrix Wtheta = a_theta * S_theta * V * D_theta1 * R_r;
  SMatrix Wphi = a_phi * S_phi * V * D_phi1 * R_r;
  SMatrix Wr = a_r * S_r * V * D_r1 * R_r;
  Ws.makeCompressed();
  Wtheta.makeCompressed();
  Wphi.makeCompressed();
  Wr.makeCompressed();

  VectorXd Ws_m0 = Ws * m0;
  // VectorXd Wtheta_m0 = Wtheta * m0;
  // VectorXd Wphi_m0 = Wphi * m0;
  // VectorXd Wr_m0 = Wr * m0;

  int nrow =
      (use_cross_gradient_constraint == true) ? (Nd + 7 * Nm) : (Nd + 4 * Nm);
  SMatrix A(nrow, Nm);
  A.reserve(Nd * Nm + Ws.nonZeros() + Wtheta.nonZeros() + Wphi.nonZeros() +
            Wr.nonZeros() + 6 * Nm);
  VectorXd b(nrow);
  A.setZero();
  b.setZero();
  // A.topRows(Nd) = (Wd * G).sparseView();
  b.head(Nd) = Wd * dobs;
  // cout<<"x"<<Nm<<endl;

  SMatrix P(Nm, Nm);
  P.reserve(Nm);
  P.setZero();

  double lambda = this->max_lambda;

  VectorXd delta_x(Nm), m_trial(Nm);  // m_trial_last(Nm), m_trial_last2(Nm);
  double misfit, misfit_last1, misfit_last2, misfit_last_iteration;

  m = m_ini;
  if (use_petrophysical_constraint) {
    m = m0;
    for (int i = 0; i < Nm; i++) {
      if (m(i) > m_max(i) || abs(m(i) - m_max(i)) < 1e-7) {
        m(i) = m_max(i) - 0.01;
      }
      if (m(i) < m_min(i) || abs(m(i) - m_min(i)) < 1e-7) {
        m(i) = m_min(i) + 0.01;
      }
    }
  }

  // cout<<G.rows()<<","<<G.cols()<<endl;
  // cout<<Wd.rows()<<","<<Wd.cols()<<endl;
  // cout<<m.rows()<<endl;
  // cout<<dobs.rows()<<endl;
  // cout<<dobs<<endl;
  misfit = sqrt((Wd * (G * m - dobs)).squaredNorm() / Nd);

  misfit_last_iteration = 5 * misfit;

  for (int id = 0; id < Nm; id++) {
    P.coeffRef(id, id) =
        (m_max(id) - m(id)) * (m(id) - m_min(id)) / (m_max(id) - m_min(id));
  }
  LeastSquaresConjugateGradient<SMatrix> lscg;
  lscg.setTolerance(cg_tol);

  ofstream out_lambda_misfit("lambda_misfit_GN");
  out_lambda_misfit << setw(25) << "lambda" << setw(25) << "misfit" << endl;

  ofstream out_iteration_misfit("Iteration_misfit_GN");
  out_iteration_misfit << setw(25) << "Iteranation number" << setw(25)
                       << "number of tried regulartion parameters" << setw(25)
                       << "misfit" << endl;

  int i, j;
  double lambda_opt;
  double misfit_opt;
  VectorXd m_opt;
  misfit_opt = misfit;
  m_opt = m;

  int c_refinement = 0;
  int i_between_refinement = 0;
  int nnz = (use_cross_gradient_constraint == true)
                ? (Nd * Nm + Ws.nonZeros() + Wtheta.nonZeros() +
                   Wphi.nonZeros() + Wr.nonZeros() + T_r.nonZeros() +
                   T_theta.nonZeros() + T_phi.nonZeros())
                : (Nd * Nm + Ws.nonZeros() + Wtheta.nonZeros() +
                   Wphi.nonZeros() + Wr.nonZeros());
  A.reserve(nnz);
  cout << "Set G to A" << endl;
  // A.topRows(Nd) = (Wd * G).sparseView();
  b.head(Nd) = Wd * dobs;

  for (size_t r_id = 0; r_id < Nd; r_id++) {
    for (size_t c_id = 0; c_id < Nm; c_id++) {
      A.insert(r_id, c_id) = Wd.coeffRef(r_id, r_id) * G(r_id, c_id);
    }
  }
  for (i = 0; i < GN_iter; i++) {
    cout << "The " << i + 1 << " th"
         << " Gauss-Newton Iteration:" << endl;
    cout << "Number of elements: " << Nm << endl;
    m_trial = m;
    misfit_last1 = 10 * misfit;   // misfit for last lambda
    misfit_last2 = 100 * misfit;  // misfit for the lambda before last lambda
    misfit_last_iteration = misfit;
    lambda = this->max_lambda;

    A.middleRows(Nd, Nm) = sqrt(lambda) * Ws;
    A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * Wr;
    A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wtheta;
    A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wphi;
    // cout << S_crg.rows() << endl;
    // cout << S_crg.cols() << endl;
    if (use_cross_gradient_constraint == true) {
      A.middleRows(Nd + 4 * Nm, Nm) = this->a_crg * S_crg * T_r;
      A.middleRows(Nd + 5 * Nm, Nm) = this->a_crg * S_crg * T_theta;
      A.middleRows(Nd + 6 * Nm, Nm) = this->a_crg * S_crg * T_phi;
    }
    A.makeCompressed();

    b.segment(Nd, Nm) = sqrt(lambda) * Ws_m0;
    // b.segment(Nd + Nm, Nm) = sqrt(lambda) * Wr_m0;
    // b.segment(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wtheta_m0;
    // b.segment(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wphi_m0;

    if (cg_iteration_factor < 1) {
      lscg.setMaxIterations(cg_iteration_factor * Nm);
    } else {
      lscg.setMaxIterations(int(cg_iteration_factor));
    }
    for (int id = 0; id < Nm; id++) {
      int n = 1.0;
      P.coeffRef(id, id) = n * (m_max(id) - m(id)) * (m(id) - m_min(id)) /
                           (m_max(id) - m_min(id));
    }
    delta_x.setZero();

    vector<double> curvature, phim_vec, phid_vec, lambda_series;
    int n_count = 0;
    for (j = 0; j < n_lambda; j++) {
      lscg.compute(A * P);
      // delta_x = lscg.solveWithGuess(b - A * m, delta_x);
      delta_x = lscg.solve(b - A * m);

      // m_trial_last2 = m_trial_last;
      // m_trial_last = m_trial;
      // double temp, temp2, temp3;
      // update
      // #pragma omp parallel for
      for (int id = 0; id < Nm; id++) {
        double temp = exp(delta_x(id));
        double temp2 = (m_min(id) * (m_max(id) - m(id)) +
                        m_max(id) * (m(id) - m_min(id)) * temp);
        double temp3 = ((m_max(id) - m(id)) + (m(id) - m_min(id)) * temp);
        if (std::isinf(temp) || std::isinf(temp2)) {
          m_trial(id) = m_max(id) - 1e-9;
        } else if (fabs(temp3) < 1e-10) {
          m_trial(id) = m_min(id) + 1e-9;
        } else {
          m_trial(id) = temp2 / temp3;
        }
      }
      n_count++;
      double phid = ( Wd*(G * m_trial - dobs)).squaredNorm();
      double phim =
          (Ws * m_trial - m0).squaredNorm() + (Wr * m_trial).squaredNorm() +
          (Wtheta * m_trial).squaredNorm() + (Wphi * m_trial).squaredNorm();
      // double phim=m_trial.squaredNorm();
      phid_vec.push_back(phid);
      phim_vec.push_back(phim);
      lambda_series.push_back(lambda);

      misfit_last2 = misfit_last1;
      misfit_last1 = misfit;
      misfit = sqrt((Wd * (G * m_trial - dobs)).squaredNorm() / Nd);
      // cout << "  Lambda=" << lambda << ", ";
      std::cout << "  #" << j << " lambda=" << lambda << ", misfit=" << misfit
                << " (CG iterations: " << lscg.iterations() << ","
                << " error: " << lscg.error() << ")" << std::endl;
      out_lambda_misfit << scientific << setw(25) << setprecision(14) << lambda
                        << scientific << setw(25) << setprecision(14)
                        << (misfit) << endl;

      if (j == 0 || misfit < misfit_opt) {
        misfit_opt = misfit;
        m_opt = m_trial;
        lambda_opt = lambda;
      }
      if ((std::abs(misfit - misfit_last1) / min(misfit, misfit_last1) <
           stag_tol) &&
          (std::abs(misfit_last1 - misfit_last2) /
               min(misfit_last1, misfit_last2) <
           stag_tol) &&
          j > 1) {
        cout << "  Misfit stagnates, go to next Gauss-Newton iteration."
             << endl;
        break;
      }

      if (misfit > misfit_last1 && misfit_last1 > misfit_last2 && j > 1) {
        cout
            << "  Misfit starts to increase, go to next Gauss-Newton iteration."
            << endl;
        break;
      } else {
        lambda = lambda * lambda_decreasing_rate;

        A.middleRows(Nd, Nm) = sqrt(lambda) * Ws;

        A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * Wr;

        A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wtheta;

        A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wphi;

        // if (use_cross_gradient_constraint == true) {
        //   A.middleRows(Nd + 4 * Nm, Nm) = this->a_crg * S_crg * T_r;
        //   A.middleRows(Nd + 5 * Nm, Nm) = this->a_crg * S_crg * T_theta;
        //   A.middleRows(Nd + 6 * Nm, Nm) = this->a_crg * S_crg * T_phi;
        // }

        b.segment(Nd, Nm) = sqrt(lambda) * Ws_m0;
        // b.segment(Nd + Nm, Nm) = sqrt(lambda) * Wr_m0;
        // b.segment(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wtheta_m0;
        // b.segment(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wphi_m0;
      }
      // if (misfit < target_misfit) {
      //   break;
      // }
    }
    // lambda=lambda_opt;
    // m = m_opt;
    // misfit = misfit_opt;

    // calculating curvature
    double max_phi_m = phim_vec[0], max_phi_d = phid_vec[0];
    for (int i_temp = 0; i_temp < n_count; i_temp++) {
      if (phim_vec[i_temp] > max_phi_m)
        max_phi_m = phim_vec[i_temp];
      if (phid_vec[i_temp] > max_phi_d)
        max_phi_d = phid_vec[i_temp];
    }
    for (int i_temp = 0; i_temp < n_count; i_temp++) {
      phim_vec[i_temp] = phim_vec[i_temp] / max_phi_m;
      phid_vec[i_temp] = phid_vec[i_temp] / max_phi_d;
    }
    // vector<double> log_phim_vec, log_phid_vec;
    // for (int i_temp = 0; i_temp < n_count; i_temp++) {
    //   log_phim_vec.push_back(log(phim_vec[i_temp])) ;
    //   log_phid_vec.push_back(log(phid_vec[i_temp])) ;
    // }    
    get_curvature(phim_vec, phid_vec, curvature);
    // get_curvature(log_phim_vec, log_phid_vec, curvature);

    ofstream Lcurve(string("L_curve") + to_string(i));
    for (int i_temp = 0; i_temp < n_count; i_temp++) {
      Lcurve << scientific << setprecision(15) << setw(25) << left
             << lambda_series[i_temp];
      Lcurve << scientific << setprecision(15) << setw(25) << left
             << log(lambda_series[i_temp]);
      Lcurve << scientific << setprecision(15) << setw(25) << left
             << phid_vec[i_temp];
      Lcurve << scientific << setprecision(15) << setw(25) << left
             << phim_vec[i_temp];
      Lcurve << scientific << setprecision(15) << setw(25) << left
             << curvature[i_temp] << endl;
    }
    unsigned int max_index = 0;
    for (int i_temp = 1; i_temp < curvature.size(); i_temp++) {
      if (curvature[i_temp] > curvature[max_index]) {
        max_index = i_temp;
      }
    }
    double max_curv_lamb = lambda_series[max_index];
    ofstream Lcurve_max_point(string("L_curve_max_point") + to_string(i));
    Lcurve_max_point << scientific << setprecision(15) << setw(25) << left
                     << lambda_series[max_index];
    Lcurve_max_point << scientific << setprecision(15) << setw(25) << left
                     << log(lambda_series[max_index]);
    Lcurve_max_point << scientific << setprecision(15) << setw(25) << left
                     << phid_vec[max_index];
    Lcurve_max_point << scientific << setprecision(15) << setw(25) << left
                     << phim_vec[max_index];
    Lcurve_max_point << scientific << setprecision(15) << setw(25) << left
                     << curvature[max_index] << endl;

    out_iteration_misfit << fixed << setw(25) << i + 1 << fixed << setw(25)
                         << ((j == n_lambda) ? n_lambda : (j + 1)) << scientific
                         << setw(25) << setprecision(14) << misfit << endl;
    // solution with optimal lambda
    lambda = max_curv_lamb;
    A.middleRows(Nd, Nm) = sqrt(lambda) * Ws;
    A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * Wr;
    A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wtheta;
    A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wphi;
    b.segment(Nd, Nm) = sqrt(lambda) * Ws_m0;
    lscg.compute(A * P);
    // delta_x = lscg.solveWithGuess(b - A * m, delta_x);
    delta_x = lscg.solve(b - A * m);
    // #pragma omp parallel for
    for (int id = 0; id < Nm; id++) {
      double temp = exp(delta_x(id));
      double temp2 = (m_min(id) * (m_max(id) - m(id)) +
                      m_max(id) * (m(id) - m_min(id)) * temp);
      double temp3 = ((m_max(id) - m(id)) + (m(id) - m_min(id)) * temp);
      if (std::isinf(temp) || std::isinf(temp2)) {
        m_trial(id) = m_max(id) - 1e-9;
      } else if (fabs(temp3) < 1e-10) {
        m_trial(id) = m_min(id) + 1e-9;
      } else {
        m_trial(id) = temp2 / temp3;
      }
    }
    misfit = sqrt((Wd * (G * m_trial - dobs)).squaredNorm() / Nd);
    m=m_trial;
    
    // cout << "  Lambda=" << lambda << ", ";
    std::cout << "  # best:" << " lambda=" << lambda << ", misfit=" << misfit
              << " (CG iterations: " << lscg.iterations() << ","
              << " error: " << lscg.error() << ")" << std::endl;

    if (record_process) {
      this->set_density_to_mesh();
      this->set_reference_model_to_mesh();

      if (use_cross_gradient_constraint == true) {
        this->mesh.out_model_netcdf(
            string("Structural_constraint_at_") + to_string(i) + string(".nc"),
            2, "crg", "");
      }
      if (use_petrophysical_constraint == true) {
        this->mesh.out_model_netcdf(string("Converted_density_model_at_") +
                                        to_string(i) + string(".nc"),
                                    1, "ref", "");
      }

      this->mesh.out_model_netcdf(string("result_at_") + to_string(i) +
                                  string(".nc"));
      this->mesh.out_model_vtk(string("result_at_") + to_string(i) +
                               string(".vtk"));
    }

    if (misfit < target_misfit || i == GN_iter - 1) {
      if (misfit < target_misfit) {
        cout << "Stop iteration, because the target misift has been achieved"
             << endl;
      } else if (i == GN_iter - 1) {
        cout << "Stop iteration, because the maximum iteration has been reached"
             << endl;
      }
      cout << "Gauss-Newton iteration number:" << i + 1 << endl;
      break;
    }

    bool flag1 = abs(misfit - misfit_last_iteration) /
                     min(misfit, misfit_last_iteration) <
                 stag_tol;
    if (c_refinement >= max_refinement_number && flag1) {
      cout << "Stop iteration, because the misfit stagnated" << endl;
      cout << "Gauss-Newton iteration number:" << i + 1 << endl;
      break;
    }

    if ((c_refinement < max_refinement_number) &&
        (++i_between_refinement % interval_between_refinements == 0 || flag1)) {
      cout << "\n";
      i_between_refinement = 0;
      ++c_refinement;
      if (flag1) {
        cout << "Misfit stagnates, refine mesh ..." << endl;
      } else {
        cout << interval_between_refinements
             << ((i_between_refinement == 1) ? string(" iteration has")
                                             : string(" iterations have"))
             << " passed since last refinement, refine mesh ..." << endl;
      }
      cout << c_refinement
           << (c_refinement == 1
                   ? ("st")
                   : (c_refinement == 2
                          ? ("nd")
                          : (c_refinement == 3 ? ("rd") : ("th"))))
           << " refinement." << endl
           << endl;

      if (use_cross_gradient_constraint && use_petrophysical_constraint) {
        refine_mesh(refinement_percentage, *interpolator_m0s, *interpolator_m0,
                    reference_surface);
        // cout<<"stop here"<<endl;
        // abort();
      } else if (use_cross_gradient_constraint &&
                 (!use_petrophysical_constraint)) {
        refine_mesh(refinement_percentage, *interpolator_m0s, reference_surface,
                    "crg");
        // cout<<"2"<<endl;
      } else if (use_petrophysical_constraint &&
                 (!use_cross_gradient_constraint)) {
        refine_mesh(refinement_percentage, *interpolator_m0, reference_surface,
                    "pet");
        // cout<<"2"<<endl;
      } else {
        refine_mesh(refinement_percentage);
        // cout<<"2"<<endl;
      }

      cout << "Finished refinement" << endl;

      // #pragma omp sections
      //       {
      // #pragma omp section
      //         {
      cout << "Initialize Ws" << endl;
      Ws.resize(Nm, Nm);
      Ws.reserve(2 * Nm);
      Ws = a_s * S_s * V * D_s * R_r;
      Ws.makeCompressed();
      Ws.data().squeeze();
      // }

      // #pragma omp section
      //         {
      cout << "Initialize Wtheta" << endl;
      Wtheta.resize(Nm, Nm);
      Wtheta.reserve(2 * Nm);
      Wtheta = a_theta * S_theta * V * D_theta1 * R_r;
      Wtheta.makeCompressed();
      Wtheta.data().squeeze();
      // }

      // #pragma omp section
      //         {
      cout << "Initialize Wphi" << endl;
      Wphi.resize(Nm, Nm);
      Wphi.reserve(2 * Nm);
      Wphi = a_phi * S_phi * V * D_phi1 * R_r;
      Wphi.makeCompressed();
      Wphi.data().squeeze();
      // }

      // #pragma omp section
      //         {
      cout << "Initialize Wr" << endl;
      Wr.resize(Nm, Nm);
      Wr.reserve(2 * Nm);
      Wr = a_r * S_r * V * D_r1 * R_r;
      Wr.makeCompressed();
      Wr.data().squeeze();
      //   }
      // }

      // update_S_crg();
      // cout<<S_crg<<endl;

      Ws_m0 = Ws * m0;
      // Wtheta_m0 = Wtheta * m0;
      // Wphi_m0 = Wphi * m0;
      // Wr_m0 = Wr * m0;
      if (use_cross_gradient_constraint == true) {
        this->set_Tmatrix();
      }

      nrow = (use_cross_gradient_constraint == true) ? (Nd + 7 * Nm)
                                                     : (Nd + 4 * Nm);
      cout << "Resizing A" << endl;
      A.resize(nrow, Nm);
      A.data().squeeze();
      VectorXi vecxi(nrow);
      for (int kk = 0; kk < Nd; kk++) {
        vecxi(kk) = Nm;
      }

      for (int kk = 0; kk < Ws.outerSize(); kk++) {
        int count_inner_nonzeros = 0;
        for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Ws, kk);
             it; ++it) {
          count_inner_nonzeros++;
        }
        vecxi(Nd + kk) = count_inner_nonzeros;
      }
      for (int kk = 0; kk < Wr.outerSize(); kk++) {
        int count_inner_nonzeros = 0;
        for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Wr, kk);
             it; ++it) {
          count_inner_nonzeros++;
        }
        // cout<<count_inner_nonzeros<<endl;
        vecxi(Nd + Nm + kk) = count_inner_nonzeros;
      }
      for (int kk = 0; kk < Wtheta.outerSize(); kk++) {
        int count_inner_nonzeros = 0;
        for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Wtheta,
                                                                     kk);
             it; ++it) {
          count_inner_nonzeros++;
        }
        vecxi(Nd + 2 * Nm + kk) = count_inner_nonzeros;
      }
      for (int kk = 0; kk < Wphi.outerSize(); kk++) {
        int count_inner_nonzeros = 0;
        for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Wphi, kk);
             it; ++it) {
          count_inner_nonzeros++;
        }
        vecxi(Nd + 3 * Nm + kk) = count_inner_nonzeros;
      }
      if (use_cross_gradient_constraint) {
        for (int kk = 0; kk < T_r.outerSize(); kk++) {
          int count_inner_nonzeros = 0;
          for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(T_r, kk);
               it; ++it) {
            count_inner_nonzeros++;
          }
          vecxi(Nd + 4 * Nm + kk) = count_inner_nonzeros;
        }
        for (int kk = 0; kk < T_theta.outerSize(); kk++) {
          int count_inner_nonzeros = 0;
          for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(T_theta,
                                                                       kk);
               it; ++it) {
            count_inner_nonzeros++;
          }
          vecxi(Nd + 5 * Nm + kk) = count_inner_nonzeros;
        }
        for (int kk = 0; kk < T_phi.outerSize(); kk++) {
          int count_inner_nonzeros = 0;
          for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(T_phi,
                                                                       kk);
               it; ++it) {
            count_inner_nonzeros++;
          }
          vecxi(Nd + 6 * Nm + kk) = count_inner_nonzeros;
        }
      }
      A.reserve(vecxi);
      cout << "Set G to A" << endl;

      for (size_t r_id = 0; r_id < Nd; r_id++) {
        for (size_t c_id = 0; c_id < Nm; c_id++) {
          A.insert(r_id, c_id) = Wd.coeffRef(r_id, r_id) * G(r_id, c_id);
        }
      }

      cout << "Resizing b" << endl;
      b.resize(nrow, 1);
      b.setZero();
      b.head(Nd) = Wd * dobs;

      P.resize(Nm, Nm);
      P.reserve(Nm);
      P.data().squeeze();

      // P.setZero();

      delta_x.resize(Nm, 1);
      m_trial.resize(Nm, 1);
    }
    // cout << i << endl;
  }
  if (i == GN_iter) {
    cout << "Gauss-Newton iteration number:" << GN_iter << endl;
  }

  this->set_density_to_mesh();

  final_lambda = lambda;
  final_misfit = misfit;
  cout << "Lambda=" << lambda << ", ";
  cout << "Misfit=" << misfit << endl;
}

void GN_Lcurve::set_CG_parameter(double cg_tol_, double cg_iteration_factor_) {
  cg_tol = cg_tol_;
  cg_iteration_factor_ = abs(cg_iteration_factor_);
  this->cg_iteration_factor = cg_iteration_factor_;
}

void GN_Lcurve::display_inversion_parameters() const {
  display_info_fields();
  cout << endl
       << GLQ_order << "th order "
       << "Gauss-Legendre quadrature is used in forward modeling " << endl;

  cout << "alpha_s=" << a_s << endl;
  cout << "alpha_r=" << a_r << endl;
  cout << "alpha_theta=" << a_theta << endl;
  cout << "alpha_phi=" << a_phi << endl;

  cout << "alpha_crg=" << a_crg << endl;

  cout << endl;
  cout << "Maximum regularizaton parameter is " << max_lambda << endl;
  cout << "Maximal number of regularization parameters (lambda) used for "
          "trials: "
       << n_lambda << endl
       << endl;

  cout << "Target misfit is " << target_misfit << endl;

  cout << "Stagnation factor is " << stag_tol
       << " (the inversion stagnate when the relative difference of misfits "
          "at "
          "2 consecutive iterations is smaller than this factor)"
       << endl;

  cout << endl;
  if (cg_iteration_factor < 1) {
    cout << "Maximal number of iterations of the conjugate gradient method is "
         << cg_iteration_factor * 100 << "% of the element number" << endl;
  } else {
    cout << "Maximal number of iterations of the conjugate gradient method is "
         << cg_iteration_factor << endl;
  }
  cout << "Tolerance value for the conjugate gradient method is " << cg_tol
       << endl;
  cout << endl;

  cout << "Number of Gauss-Newton iterations: " << GN_iter << endl;

  cout << "Model density contrast limits: "
       << "[" << m_min.minCoeff() << ", " << m_max.maxCoeff() << "]"
       << "kg/m3" << endl;

  if (max_refinement_number > 0) {
    cout << endl;
    cout << "The inversion mesh will be adaptively refined at every "
         << interval_between_refinements << " iteration." << endl;
    cout << "Maximal times of refinements: " << max_refinement_number << endl;
    cout << "Percentage of elements marked for refinement is "
         << refinement_percentage * 100 << "%" << endl;
  } else {
    cout << endl;
    cout << "The inversion mesh will not be refined." << endl;
  }
  cout << "Will the inversion model be recorded at each iteration?";
  if (record_process == 1) {
    cout << " Yes" << endl;
  } else {
    cout << " No" << endl;
  }
  cout << endl;
  cout << "Use cross-gradient constraint model? "
       << ((use_cross_gradient_constraint == true) ? ("Yes") : ("No")) << endl;
  cout << "Use petrophysical constraint model? "
       << ((use_petrophysical_constraint == true) ? ("Yes") : ("No")) << endl;
  cout << endl;
}

// void GN_Lcurve::get_curvature(const vector<double>& x,
//                                          const vector<double>& y,
//                                          vector<double>& curvature) {
//   curvature.resize(x.size());
//   for (int i = 1; i < x.size() - 1; i++) {
//     // https://zhuanlan.zhihu.com/p/72083902
//     double x1 = x[i - 1], x2 = x[i], x3 = x[i + 1];
//     double y1 = y[i - 1], y2 = y[i], y3 = y[i + 1];
//     double ta = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
//     double tb = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
//     Vector3d xv;
//     Matrix3d M;
//     xv << x1, x2, x3;
//     M << 1., -ta * ta, ta * ta, 1, 0, 0, 1, tb, tb * tb;
//     Vector3d yv;
//     yv << y1, y2, y3;
//     Vector3d a = M.colPivHouseholderQr().solve(xv);
//     Vector3d b = M.colPivHouseholderQr().solve(yv);

//     double xpp = 2 * a(2);
//     double ypp = 2 * b(2);
//     double xp = a(1);
//     double yp = b(1);
//     curvature[i] = (xpp * yp - xp * ypp) / pow(xp * xp + yp * yp, 1.5);
//     if (i == 1) {
//       xpp = 2 * a(2);
//       ypp = 2 * b(2);
//       xp = a(1) + 2 * a(2) * ta;
//       yp = b(1) + 2 * b(2) * ta;
//       curvature[0] = (xpp * yp - xp * ypp) / pow(xp * xp + yp * yp, 1.5);
//     }
//     if (i == x.size() - 2) {
//       xpp = 2 * a(2);
//       ypp = 2 * b(2);
//       xp = a(1) + 2 * a(2) * tb;
//       yp = b(1) + 2 * b(2) * tb;
//       curvature[x.size() - 1] =
//           (xpp * yp - xp * ypp) / pow(xp * xp + yp * yp, 1.5);
//     }
//   }
// }
void GN_Lcurve::get_curvature(const vector<double>& x,
                              const vector<double>& y,
                              vector<double>& curvature) {
  int npoint = x.size();
  curvature.resize(npoint);
  MatrixXd B(npoint, npoint);
  B.setZero();
  VectorXd h(npoint - 1);
  for (int i = 0; i < npoint - 1; i++) {
    h(i) = x[i + 1] - x[i];
  }
  VectorXd alpha(npoint);
  VectorXd beta(npoint);
  alpha(0) = 0;
  beta(npoint - 1) = 0;

  VectorXd d(npoint);
  d(0) = 0;
  d(npoint - 1) = 0;

  for (int i = 1; i < npoint - 1; i++) {
    double len = h(i - 1) + h(i);
    alpha(i) = h(i) / len;
    beta(i) = h(i - 1) / len;
    d(i) =
        6.0 * ((y[i + 1] - y[i]) / h(i) - (y[i] - y[i - 1]) / h(i - 1)) / len;
  }
  B(0, 0) = 2.0;
  B(0, 1) = alpha(0);
  B(npoint - 1, npoint - 2) = beta(npoint - 1);
  B(npoint - 1, npoint - 1) = 2.0;

  for (int i = 1; i < npoint - 1; i++) {
    B(i, i - 1) = beta(i);
    B(i, i) = 2.0;
    B(i, i + 1) = alpha(i);
  }
  VectorXd ypp = B.colPivHouseholderQr().solve(d);
  VectorXd yp(npoint);
  for (int i = 0; i < npoint - 1; i++) {
    yp(i) = -ypp(i) * h(i) / 3.0 - ypp(i + 1) * h(i) / 6.0 +
            (y[i + 1] - y[i]) / h(i);
  }
  yp(npoint - 1) = h(npoint - 2) * ypp(npoint - 1) / 3.0 +
                   h(npoint - 2) * ypp(npoint - 2) / 6.0 +
                   (y[npoint - 1] - y[npoint - 2]) / h(npoint - 2);

  for (int i = 0; i < npoint; i++) {
    curvature[i] =
        abs(ypp(i)) /
        sqrt((1 + yp(i) * yp(i)) * (1 + yp(i) * yp(i)) * (1 + yp(i) * yp(i)));
  }
}
