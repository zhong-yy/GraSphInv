#ifndef Crust1Correction_H
#define Crust1Correction_H
#include <fstream>
#include<iostream>
#include <string>
#include"timer.h"
#include "Fwd.h"
#include "Mesh.h"
#include "Observation.h"
class Crust1Correction {
 public:
  Crust1Correction();
  void show_crust_info();
  void read_crust_data();
  void build_mesh_of_sediments();
  void build_mesh_of_crystalline_crust();
  void build_mesh_of_moho();
  void set_ob(const Observation& ob_) { ob = ob_; }
  void gravity_effects_of_sediments();
  void gravity_effects_of_crystalline_crusts();
  void gravity_effects_of_moho();

  void line_process(std::string& line, const std::string comment_str = "#");
  istream& next_valid_line(istream& is, string& str);

  void read_data_to_be_corrected();

  void read_configuration_file(string file_name);

  void out_file();

 private:
  string path_to_crust1_0;
  vector<vector<vector<double> > >crust1_bnds;
  vector<vector<vector<double> > >crust1_rho;
  // double crust1_bnds[180][360][9];
  // double crust1_rho[180][360][9];
  double reference_surface;
  double reference_crust_density;
  double reference_mantle_density;
  // double reference_upper_crystalline_crust;
  // double reference_middle_crystalline_crust;
  // double reference_lower_crystalline_crust;
  Mesh sediments;
  Mesh crystalline_crust;
  Mesh moho;
  int n_fields;
  unsigned long long field_flag;
  Observation ob;

  VectorXd gra_sediments;
  VectorXd gra_crystalline_crust;
  VectorXd gra_moho;


  VectorXd dobs;

  string data_file;
  double height;
};
#endif