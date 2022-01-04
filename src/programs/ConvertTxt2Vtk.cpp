#include"Mesh.h"
#include<iostream>
using namespace std;
int main(int argc, char** argv)
{
  if (argc < 3) {
    printf("Usage: %s input_txt_file_name output_vtk_filename\n", argv[0]);
    return 1;
  }
  string txtfile(argv[1]);
  string vtkfile(argv[2]);

  Mesh mesh;
  mesh.convert_txt_2_vtk(txtfile,vtkfile);
  return 0;
}
