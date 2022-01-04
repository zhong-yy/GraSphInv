#include <cassert>
#include <iostream>
#include <string>
#include "Crust1Correction2.h"
using namespace std;
int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Usage: %s input_paramters_filename\n", argv[0]);
    return 1;
  }
  assert(argc == 2);
  string configuration_file(argv[1]);
  Crust1Correction2 crust_correction;

  crust_correction.read_configuration_file(configuration_file);
  crust_correction.read_data_to_be_corrected();
  cout << "A" << endl;
  crust_correction.read_crust_data();
  cout << "B" << endl;
  crust_correction.build_mesh_of_sediments();
  cout << "C" << endl;
  crust_correction.build_mesh_of_crystalline_crust();
  crust_correction.build_mesh_of_moho();
//   cout << "D" << endl;
  crust_correction.gravity_effects_of_sediments();
  crust_correction.gravity_effects_of_crystalline_crusts();
  crust_correction.gravity_effects_of_moho();
  crust_correction.out_file();
  return 0;
}
