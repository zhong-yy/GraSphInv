#ifndef Crust1Correction2_H
#define Crust1Correction2_H
#include <fstream>
#include<iostream>
#include <string>
#include"timer.h"
#include "Fwd.h"
#include "Mesh.h"
#include "Observation.h"
class Crust1Correction2 {
 public:
  Crust1Correction2();
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
  vector<vector<vector<double> > >crust1_bnds;//boundaries
  vector<vector<vector<double> > >crust1_rho;
  // double crust1_bnds[180][360][9];
  // double crust1_rho[180][360][9];
  double reference_surface;
  double reference_upper_crust;
  double reference_lower_crust;
  double reference_upper_crust_depth;
  double reference_lower_crust_depth;

  double reference_mantle;
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