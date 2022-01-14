#include "Crust1Correction.h"
Crust1Correction::Crust1Correction() {
    reference_surface = 6378137;
    reference_upper_crust = 2700;
    reference_lower_crust = 2940;
    reference_upper_crust_depth = 15000;  // 15km
    reference_lower_crust_depth = 40000;  // 40km
    // reference_mantle = 3357;
    reference_mantle = 3200;
    path_to_crust1_0 = "";
    n_fields = 0;
    field_flag = 0;
    crust1_bnds.resize(180);
    crust1_rho.resize(180);
    for (int i = 0; i < 180; i++) {
        crust1_bnds[i].resize(360);
        crust1_rho[i].resize(360);
        for (int j = 0; j < 360; j++) {
            crust1_bnds[i][j].resize(9);
            crust1_rho[i][j].resize(9);
        }
    }
}
void Crust1Correction::read_crust_data() {
    string filename_bnds = path_to_crust1_0 + "/crust1.bnds";
    string filename_rho = path_to_crust1_0 + "/crust1.rho";
    cout << "Reading " << filename_bnds << endl;
    cout << "Reading " << filename_rho << endl;
    ifstream ifs1(filename_bnds.c_str());
    ifstream ifs2(filename_rho.c_str());
    assert(ifs1.good());
    assert(ifs2.good());
    // first point's coordinate is (179.5W,89.5N)
    for (int i = 0; i < 180; i++) {
        for (int j = 0; j < 360; j++) {
            for (int k = 0; k < 9; k++) {
                ifs1 >> crust1_bnds[i][j][k];
                ifs2 >> crust1_rho[i][j][k];
            }
        }
    }

    show_crust_info();
}

void Crust1Correction::show_crust_info() {
    double average_upper_crystalline = 0;
    double average_middle_crystalline = 0;
    double average_lower_crystalline = 0;
    int num_upper_crust = 0;
    int num_middle_crust = 0;
    int num_lower_crust = 0;
    for (int i = 0; i < 180; i++) {
        for (int j = 0; j < 360; j++) {
            // upper crust
            if (abs(crust1_bnds[i][j][5] - crust1_bnds[i][j][6]) > 1e-5) {
                average_upper_crystalline +=
                    crust1_rho[i][j][5] * 1000;  // convert g/cm3 to kg/m3
                num_upper_crust++;
            }
            // middle crust
            if (abs(crust1_bnds[i][j][6] - crust1_bnds[i][j][7]) > 1e-5) {
                average_middle_crystalline += crust1_rho[i][j][6] * 1000;
                num_middle_crust++;
            }
            if (abs(crust1_bnds[i][j][7] - crust1_bnds[i][j][8]) > 1e-5) {
                average_lower_crystalline += crust1_rho[i][j][7] * 1000;
                num_lower_crust++;
            }
        }
    }
    cout << "num_upper_crust=" << num_upper_crust << endl;
    cout << "num_middle_crust=" << num_middle_crust << endl;
    cout << "num_lower_crust=" << num_lower_crust << endl;
    average_upper_crystalline = average_upper_crystalline / num_upper_crust;
    average_middle_crystalline = average_middle_crystalline / num_middle_crust;
    average_lower_crystalline = average_lower_crystalline / num_lower_crust;
    cout << "average_upper_crystalline:" << average_upper_crystalline << endl;
    cout << "average_middle_crystalline:" << average_middle_crystalline << endl;
    cout << "average_lower_crystalline:" << average_lower_crystalline << endl;
    double deepest_upper_sediment = 0;
    double deepest_middle_sediment = 0;
    double deepest_lower_sediment = 0;
    double deepest_upper_crust = 0;
    double deepest_middle_crust = 0;
    double deepest_lower_crust = 0;
    double shollowest_lower_crust = abs(crust1_bnds[0][0][8]);
    for (int i = 0; i < 180; i++) {
        for (int j = 0; j < 360; j++) {
            if (abs(crust1_bnds[i][j][3]) > deepest_upper_sediment) {
                deepest_upper_sediment = abs(crust1_bnds[i][j][3]);
            }
            if (abs(crust1_bnds[i][j][4]) > deepest_middle_sediment) {
                deepest_middle_sediment = abs(crust1_bnds[i][j][4]);
            }
            if (abs(crust1_bnds[i][j][5]) > deepest_lower_sediment) {
                deepest_lower_sediment = abs(crust1_bnds[i][j][5]);
            }
            if (abs(crust1_bnds[i][j][6]) > deepest_upper_crust) {
                deepest_upper_crust = abs(crust1_bnds[i][j][6]);
            }
            if (abs(crust1_bnds[i][j][7]) > deepest_middle_crust) {
                deepest_middle_crust = abs(crust1_bnds[i][j][7]);
            }
            if (abs(crust1_bnds[i][j][8]) > deepest_lower_crust) {
                deepest_lower_crust = abs(crust1_bnds[i][j][8]);
            }

            if (abs(crust1_bnds[i][j][8]) < shollowest_lower_crust) {
                shollowest_lower_crust = abs(crust1_bnds[i][j][8]);
            }
        }
    }

    cout << "deepest_upper_sediment=" << deepest_upper_sediment << endl;
    cout << "deepest_middle_sediment=" << deepest_middle_sediment << endl;
    cout << "deepest_lower_sediment=" << deepest_lower_sediment << endl;
    cout << "deepest_upper_crust=" << deepest_upper_crust << endl;
    cout << "deepest_middle_crust=" << deepest_middle_crust << endl;
    cout << "deepest_lower_crust=" << deepest_lower_crust << endl;
    cout << "shallowest_lower_crust=" << shollowest_lower_crust << endl;
}

void Crust1Correction::build_mesh_of_sediments() {
    sediments.cells.push_back(vector<Cell*>(0));
    int id = 0;
    double r_upper_crust = reference_surface - reference_upper_crust_depth;
    // deepest sediment is 21km
    for (int i = 0; i < 180; i++) {
        for (int j = 0; j < 360; j++) {
            // upper sediment
            double r0, theta0, phi0;
            double r1, theta1, phi1;
            theta0 = (90.0 - (90.0 - i)) * GS::PI / 180.0;
            theta1 = (90.0 - (89.0 - i)) * GS::PI / 180.0;
            phi0 = (180.0 + (-180.0 + j)) * GS::PI / 180.0;
            phi1 = (180.0 + (-179.0 + j)) * GS::PI / 180.0;
            if (abs(crust1_bnds[i][j][2] - crust1_bnds[i][j][3]) > 1e-5) {
                r0 = reference_surface + crust1_bnds[i][j][3] * 1000;
                r1 = reference_surface + crust1_bnds[i][j][2] * 1000;
                assert(r1 > r0);
                if (r0 > r_upper_crust || abs(r0 - r_upper_crust) < 1e-5) {
                    // within upper crust
                    sediments.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][2] * 1000 - reference_upper_crust);
                    id++;
                } else if (r1 < r_upper_crust ||
                           abs(r1 - r_upper_crust) < 1e-5) {
                    // within lower crust
                    sediments.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][2] * 1000 - reference_lower_crust);
                    id++;
                } else {
                    // span upper crust and lower crust
                    assert(r0 < r_upper_crust);
                    assert(r1 > r_upper_crust);
                    sediments.cells[0].push_back(new Cell(r0, theta0, phi0,
                                                          r_upper_crust, theta1,
                                                          phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][2] * 1000 - reference_lower_crust);
                    id++;

                    sediments.cells[0].push_back(new Cell(r_upper_crust, theta0,
                                                          phi0, r1, theta1,
                                                          phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][2] * 1000 - reference_upper_crust);
                    id++;
                }
            }

            // middle sediment
            if (abs(crust1_bnds[i][j][3] - crust1_bnds[i][j][4]) > 1e-5) {
                r0 = reference_surface + crust1_bnds[i][j][4] * 1000;
                r1 = reference_surface + crust1_bnds[i][j][3] * 1000;
                assert(r1 > r0);
                if (r0 > r_upper_crust || abs(r0 - r_upper_crust) < 1e-5) {
                    // within upper crust
                    sediments.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][3] * 1000 - reference_upper_crust);
                    id++;
                } else if (r1 < r_upper_crust ||
                           abs(r1 - r_upper_crust) < 1e-5) {
                    // within lower crust
                    sediments.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][3] * 1000 - reference_lower_crust);
                    id++;
                } else {
                    // span upper crust and lower crust
                    assert(r0 < r_upper_crust);
                    assert(r1 > r_upper_crust);
                    sediments.cells[0].push_back(new Cell(r0, theta0, phi0,
                                                          r_upper_crust, theta1,
                                                          phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][3] * 1000 - reference_lower_crust);
                    id++;

                    sediments.cells[0].push_back(new Cell(r_upper_crust, theta0,
                                                          phi0, r1, theta1,
                                                          phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][3] * 1000 - reference_upper_crust);
                    id++;
                }
            }

            // lower sediment
            if (abs(crust1_bnds[i][j][4] - crust1_bnds[i][j][5]) > 1e-5) {
                r0 = reference_surface + crust1_bnds[i][j][5] * 1000;
                r1 = reference_surface + crust1_bnds[i][j][4] * 1000;
                if (r0 > r_upper_crust || abs(r0 - r_upper_crust) < 1e-5) {
                    // within upper crust
                    sediments.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][4] * 1000 - reference_upper_crust);
                    id++;
                } else if (r1 < r_upper_crust ||
                           abs(r1 - r_upper_crust) < 1e-5) {
                    // within lower crust
                    sediments.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][4] * 1000 - reference_lower_crust);
                    id++;
                } else {
                    // span upper crust and lower crust
                    assert(r0 < r_upper_crust);
                    assert(r1 > r_upper_crust);
                    sediments.cells[0].push_back(new Cell(r0, theta0, phi0,
                                                          r_upper_crust, theta1,
                                                          phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][4] * 1000 - reference_lower_crust);
                    id++;

                    sediments.cells[0].push_back(new Cell(r_upper_crust, theta0,
                                                          phi0, r1, theta1,
                                                          phi1, 0, 1, true));
                    sediments.cells[0][id]->set_parameter(
                        crust1_rho[i][j][4] * 1000 - reference_upper_crust);
                    id++;
                }
            }
        }
    }
    sediments.leaf_cells.resize(sediments.cells[0].size());
    sediments.num_leaf_cells = sediments.cells[0].size();
    for (int i = 0; i < sediments.num_leaf_cells; i++) {
        sediments.leaf_cells[i] = sediments.cells[0][i];
    }
    sediments.out_model_vtk("sediments.vtk");
}

void Crust1Correction::build_mesh_of_crystalline_crust() {
    crystalline_crust.cells.push_back(vector<Cell*>(0));
    double r_upper_crust = reference_surface - reference_upper_crust_depth;
    int id = 0;
    for (int i = 0; i < 180; i++) {
        for (int j = 0; j < 360; j++) {
            // upper crust
            double r0, theta0, phi0;
            double r1, theta1, phi1;
            theta0 = (90.0 - (90.0 - i)) * GS::PI / 180.0;
            theta1 = (90.0 - (89.0 - i)) * GS::PI / 180.0;
            phi0 = (180.0 + (-180.0 + j)) * GS::PI / 180.0;
            phi1 = (180.0 + (-179.0 + j)) * GS::PI / 180.0;
            if (abs(crust1_bnds[i][j][5] - crust1_bnds[i][j][6]) > 1e-5) {
                r0 = reference_surface + crust1_bnds[i][j][6] * 1000;
                r1 = reference_surface + crust1_bnds[i][j][5] * 1000;
                if (r0 > r_upper_crust || abs(r0 - r_upper_crust) < 1e-5) {
                    // within upper crust
                    crystalline_crust.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][5] * 1000 - reference_upper_crust);
                    id++;
                } else if (r1 < r_upper_crust ||
                           abs(r1 - r_upper_crust) < 1e-5) {
                    // within lower crust
                    crystalline_crust.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][5] * 1000 - reference_lower_crust);
                    id++;
                } else {
                    // span upper crust and lower crust
                    assert(r0 < r_upper_crust);
                    assert(r1 > r_upper_crust);
                    crystalline_crust.cells[0].push_back(
                        new Cell(r0, theta0, phi0, r_upper_crust, theta1, phi1,
                                 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][5] * 1000 - reference_lower_crust);
                    id++;

                    crystalline_crust.cells[0].push_back(
                        new Cell(r_upper_crust, theta0, phi0, r1, theta1, phi1,
                                 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][5] * 1000 - reference_upper_crust);
                    id++;
                }
            }

            // middle crust
            if (abs(crust1_bnds[i][j][6] - crust1_bnds[i][j][7]) > 1e-5) {
                r0 = reference_surface + crust1_bnds[i][j][7] * 1000;
                r1 = reference_surface + crust1_bnds[i][j][6] * 1000;
                if (r0 > r_upper_crust || abs(r0 - r_upper_crust) < 1e-5) {
                    // within upper crust
                    crystalline_crust.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][6] * 1000 - reference_upper_crust);
                    id++;
                } else if (r1 < r_upper_crust ||
                           abs(r1 - r_upper_crust) < 1e-5) {
                    // within lower crust
                    crystalline_crust.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][6] * 1000 - reference_lower_crust);
                    id++;
                } else {
                    // span upper crust and lower crust
                    assert(r0 < r_upper_crust);
                    assert(r1 > r_upper_crust);
                    crystalline_crust.cells[0].push_back(
                        new Cell(r0, theta0, phi0, r_upper_crust, theta1, phi1,
                                 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][6] * 1000 - reference_lower_crust);
                    id++;

                    crystalline_crust.cells[0].push_back(
                        new Cell(r_upper_crust, theta0, phi0, r1, theta1, phi1,
                                 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][6] * 1000 - reference_upper_crust);
                    id++;
                }
            }

            // lower crust
            if (abs(crust1_bnds[i][j][7] - crust1_bnds[i][j][8]) > 1e-5) {
                r0 = reference_surface + crust1_bnds[i][j][8] * 1000;
                r1 = reference_surface + crust1_bnds[i][j][7] * 1000;
                if (r0 > r_upper_crust || abs(r0 - r_upper_crust) < 1e-5) {
                    // within upper crust
                    crystalline_crust.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][7] * 1000 - reference_upper_crust);
                    id++;
                } else if (r1 < r_upper_crust ||
                           abs(r1 - r_upper_crust) < 1e-5) {
                    // within lower crust
                    crystalline_crust.cells[0].push_back(new Cell(
                        r0, theta0, phi0, r1, theta1, phi1, 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][7] * 1000 - reference_lower_crust);
                    id++;
                } else {
                    // span upper crust and lower crust
                    assert(r0 < r_upper_crust);
                    assert(r1 > r_upper_crust);
                    crystalline_crust.cells[0].push_back(
                        new Cell(r0, theta0, phi0, r_upper_crust, theta1, phi1,
                                 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][7] * 1000 - reference_lower_crust);
                    id++;

                    crystalline_crust.cells[0].push_back(
                        new Cell(r_upper_crust, theta0, phi0, r1, theta1, phi1,
                                 0, 1, true));
                    crystalline_crust.cells[0][id]->set_parameter(
                        crust1_rho[i][j][7] * 1000 - reference_upper_crust);
                    id++;
                }
            }
        }
    }
    crystalline_crust.leaf_cells.resize(crystalline_crust.cells[0].size());
    crystalline_crust.num_leaf_cells = crystalline_crust.cells[0].size();
    for (int i = 0; i < crystalline_crust.num_leaf_cells; i++) {
        crystalline_crust.leaf_cells[i] = crystalline_crust.cells[0][i];
    }
    crystalline_crust.out_model_vtk("crystalline_crust.vtk");
}

void Crust1Correction::build_mesh_of_moho() {
    moho.cells.push_back(vector<Cell*>(0));
    int id = 0;
    double r_upper_crust = reference_surface - reference_upper_crust_depth;
    double r_lower_crust = reference_surface - reference_lower_crust_depth;

    for (int i = 0; i < 180; i++) {
        for (int j = 0; j < 360; j++) {
            // upper crust
            double r0, theta0, phi0;
            double r1, theta1, phi1;
            theta0 = (90.0 - (90.0 - i)) * GS::PI / 180.0;
            theta1 = (90.0 - (89.0 - i)) * GS::PI / 180.0;
            phi0 = (180.0 + (-180.0 + j)) * GS::PI / 180.0;
            phi1 = (180.0 + (-179.0 + j)) * GS::PI / 180.0;
            double r_moho = reference_surface + crust1_bnds[i][j][8] * 1000;
            if (abs(r_moho - r_lower_crust) < 1e-5) {
                // moho is coincide with the reference crust-mantle boundary
                // do nothing
            } else if (r_moho > r_lower_crust) {
                if (r_moho > r_upper_crust) {
                    // moho is above the reference crust-mantle boundary
                    moho.cells[0].push_back(new Cell(r_lower_crust, theta0,
                                                     phi0, r_upper_crust,
                                                     theta1, phi1, 0, 1, true));
                    moho.cells[0][id]->set_parameter(reference_mantle -
                                                     reference_lower_crust);
                    id++;

                    moho.cells[0].push_back(new Cell(r_upper_crust, theta0,
                                                     phi0, r_moho, theta1, phi1,
                                                     0, 1, true));
                    moho.cells[0][id]->set_parameter(reference_mantle -
                                                     reference_upper_crust);
                    id++;
                } else {
                    moho.cells[0].push_back(new Cell(r_lower_crust, theta0,
                                                     phi0, r_moho, theta1, phi1,
                                                     0, 1, true));
                    moho.cells[0][id]->set_parameter(reference_mantle -
                                                     reference_lower_crust);
                    id++;
                }
            } else {
                // moho is below the reference crust-mantle boundary
                moho.cells[0].push_back(new Cell(r_moho, theta0, phi0,
                                                 r_lower_crust, theta1, phi1, 0,
                                                 1, true));
                moho.cells[0][id]->set_parameter(reference_lower_crust -
                                                 reference_mantle);
                id++;
            }
        }
    }
    moho.leaf_cells.resize(moho.cells[0].size());
    moho.num_leaf_cells = moho.cells[0].size();
    for (int i = 0; i < moho.num_leaf_cells; i++) {
        moho.leaf_cells[i] = moho.cells[0][i];
    }
    moho.out_model_vtk("moho.vtk");
}
void Crust1Correction::read_data_to_be_corrected() {
    cout << "Reading " << data_file << endl;
    ifstream input_stream(data_file.c_str());
    if (!input_stream.good()) {
        cout << "File " << data_file << "does not exist" << endl;
        abort();
    }
    assert(input_stream.good());
    string line;
    int n_obs = 0;  // number of observation points
    vector<vector<double> > g_data;
    g_data.resize(n_fields);
    while (std::getline(input_stream, line)) {
        line_process(line);
        if (line.empty()) {
            continue;
        } else {
            n_obs++;
            double r0, theta0, phi0;
            std::istringstream iss(line);
            iss >> theta0 >> phi0;
            theta0 = (90.0 - theta0) * GS::PI / 180.0;
            phi0 = (180.0 + phi0) * GS::PI / 180;

            r0 = reference_surface + height;
            ob.add_point(r0, theta0, phi0);
            double gtheta, gphi, gr;
            for (int j = 0; j < n_fields; j++) {
                double tmp;
                iss >> tmp;
                g_data[j].push_back(tmp);
                // input_stream >> dobs(n_line + j * n_obs);
            }
        }
    }
    dobs.resize(n_obs * n_fields);
    for (int j = 0; j < n_fields; j++) {
        assert(g_data[j].size() == n_obs);
        for (int i = 0; i < n_obs; i++) {
            dobs(i + j * n_obs) = g_data[j][i];
        }
    }
    cout << "Number of observation points: " << n_obs << endl;
}

void Crust1Correction::out_file() {
    bitset<32> bitvec(field_flag);
    assert(n_fields == bitvec.count());
    vector<unsigned int> field_label;
    for (int i = 0; i < 10; i++) {
        if (bitvec[i]) {
            field_label.push_back(i);
        }
    }
    // for (int i = 0; i < field_label.size(); i++) {
    //   cout << field_label[i] << endl;
    // }
    vector<string> strs = {"V",          "g_r",      "g_theta", "g_phi",
                           "T_rr",       "T_rtheta", "T_rphi",  "T_thetatheta",
                           "T_thetaphi", "T_phiphi"};

    int n_com = bitvec.count();  // number of used components
    int n_ob = ob.get_n_obs();
    int nd = gra_sediments.rows();

    assert(gra_sediments.rows() == gra_crystalline_crust.rows());

    assert(nd % n_ob == 0);
    assert(nd / n_ob == n_com);

    for (int j = 0; j < n_com; j++) {
        string file_name = "sediments_" + strs[field_label[j]];
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
            out_s << setw(30) << setprecision(15) << left
                  << gra_sediments(i + j * n_ob) << endl;
        }
    }

    for (int j = 0; j < n_com; j++) {
        string file_name = "crystalline_crust_" + strs[field_label[j]];
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
            out_s << setw(30) << setprecision(15) << left
                  << gra_crystalline_crust(i + j * n_ob) << endl;
        }
    }

    for (int j = 0; j < n_com; j++) {
        string file_name = "moho_" + strs[field_label[j]];
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
            out_s << setw(30) << setprecision(15) << left
                  << gra_moho(i + j * n_ob) << endl;
        }
    }

    for (int j = 0; j < n_com; j++) {
        string file_name = "rm_sedim_" + strs[field_label[j]];
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
            out_s << setw(30) << setprecision(15) << left
                  << rm_sedim_gra(i + j * n_ob) << endl;
        }
    }
    for (int j = 0; j < n_com; j++) {
        string file_name = "rm_sedim_cryst_" + strs[field_label[j]];
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
            out_s << setw(30) << setprecision(15) << left
                  << rm_sedim_cryst_gra(i + j * n_ob) << endl;
        }
    }
    for (int j = 0; j < n_com; j++) {
        string file_name = "mantle_" + strs[field_label[j]];
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
            out_s << setw(30) << setprecision(15) << left
                  << mantle_gra(i + j * n_ob) << endl;
        }
    }
}
void Crust1Correction::read_configuration_file(string file_name) {
    ifstream ifs(file_name.c_str());
    assert(ifs.good());
    string line;
    next_valid_line(ifs, line);
    istringstream iss(line);
    iss >> path_to_crust1_0;
    cout << path_to_crust1_0 << endl;

    next_valid_line(ifs, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> data_file;

    next_valid_line(ifs, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> height;
    cout << "Height of data from the reference surface  is " << height << endl;

    next_valid_line(ifs, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> reference_surface;
    cout << "Reference level is " << setprecision(7) << reference_surface
         << endl;

    next_valid_line(ifs, line);
    iss.clear();
    iss.str("");
    iss.str(line);

    iss >> n_fields;
    cout << endl;
    cout << n_fields << " gravity component"
         << ((n_fields == 1) ? (" is") : ("s are")) << " used" << endl;

    field_flag = 0;
    vector<unsigned int> field_label(n_fields);
    next_valid_line(ifs, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    for (int i = 0; i < n_fields; i++) {
        iss >> field_label[i];
        switch (field_label[i]) {
            case 0:
                field_flag = field_flag | Compute_V;
                break;
            case 1:
                field_flag = field_flag | Compute_g_r;
                break;
            case 2:
                field_flag = field_flag | Compute_g_theta;
                break;
            case 3:
                field_flag = field_flag | Compute_g_phi;
                break;
            case 4:
                field_flag = field_flag | Compute_T_rr;
                break;
            case 5:
                field_flag = field_flag | Compute_T_rtheta;
                break;
            case 6:
                field_flag = field_flag | Compute_T_rphi;
                break;
            case 7:
                field_flag = field_flag | Compute_T_thetatheta;
                break;
            case 8:
                field_flag = field_flag | Compute_T_thetaphi;
                break;
            case 9:
                field_flag = field_flag | Compute_T_phiphi;
                break;
        }
    }
}
void Crust1Correction::line_process(std::string& line,
                                    const std::string comment_str) {
    for (char& c : line)  // C++11以上版本的语法
    {
        //制表符 tab，逗号，分号都当作有效的分隔符，统一转成空格
        //为了避免错误，回车符和换行符也转为空格（否则无法处理空行）
        if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n')
            c = ' ';
    }

    //查找注释符所在位置，如果不存在，则得到string::npos
    int n_comment_start = line.find_first_of(comment_str);
    if (n_comment_start != std::string::npos)  //这一句必须的
        line.erase(n_comment_start);           //删除注释

    line.erase(0, line.find_first_not_of(" "));  //删除行首空格
    line.erase(line.find_last_not_of(" ") + 1);  //删除行末空格
    //调用的string& erase (size_t pos = 0, size_t len = npos);
    // len为默认参数,
    // size_t可以当作无符号整数，npos是string内置的静态常量，为size_t的最大值

    /****************************************************************
     *处理完毕。如果这一行只有空格，制表符 tab，注释，那么处理后line为空；
     *如果行首有多个空格(或者空格和tab交错)，行尾为注释，如
     *“   a b c#坐标”
     *那么处理后字符串line的行首多个空格(和tab)和行尾注释被删掉，得到
     *“a b c”
     ****************************************************************/
}

istream& Crust1Correction::next_valid_line(istream& input_stream,
                                           string& line) {
    while (std::getline(input_stream, line)) {
        line_process(line);
        if (line.empty()) {
            continue;
        } else {
            break;
        }
    }
    return input_stream;
}

void Crust1Correction::gravity_effects_of_sediments() {
    cout << "Calculating gravity effects of density contrasts in sedimentary "
            "layers"
         << endl;
    Timer t;
    t.start();
    Fwd fwd(sediments, ob, field_flag, 8);
    // fwd.compute_G();
    // const MatrixXd& G = fwd.get_G();
    VectorXd m;
    sediments.get_model_parameter_from_mesh(m);
    gra_sediments = fwd.compute_gobs_without_G(m);
    // gra_sediments = G * m;
    t.stop();
    cout << "Finished" << endl;
    cout << "Time: " << t.getElapsedTimeInSec() << " s" << endl;
}

void Crust1Correction::gravity_effects_of_crystalline_crusts() {
    Timer t;
    t.start();
    cout << "Calculating gravity effects of density contrasts in crystalline "
            "crusts"
         << endl;
    Fwd fwd(crystalline_crust, ob, field_flag, 8);
    // fwd.compute_G();
    // const MatrixXd& G = fwd.get_G();
    VectorXd m;
    crystalline_crust.get_model_parameter_from_mesh(m);
    gra_crystalline_crust = fwd.compute_gobs_without_G(m);
    // gra_crystalline_crust = G * m;
    cout << "Finished" << endl;
    t.stop();
    cout << "Time: " << t.getElapsedTimeInSec() << " s" << endl;
}

void Crust1Correction::gravity_effects_of_moho() {
    Timer t;
    t.start();
    cout << "Calculating gravity effects of moho" << endl;
    Fwd fwd(moho, ob, field_flag, 8);
    // fwd.compute_G();
    // const MatrixXd& G = fwd.get_G();
    VectorXd m;
    moho.get_model_parameter_from_mesh(m);
    gra_moho = fwd.compute_gobs_without_G(m);
    // gra_moho = G * m;
    cout << "Finished" << endl;
    t.stop();
    cout << "Time: " << t.getElapsedTimeInSec() << " s" << endl;
}
void Crust1Correction::remove_crustal_effects() {
    this->rm_sedim_gra = dobs - gra_sediments;
    this->rm_sedim_cryst_gra = dobs - gra_sediments - gra_crystalline_crust;
    this->mantle_gra = dobs - gra_sediments - gra_crystalline_crust - gra_moho;
}
