#include <iostream>

#include "Mesh.h"
using namespace std;
int main(int argc, char** argv) {
    if (argc < 3) {
        printf("Usage: %s input_txt_file_name output_vtk_filename\n", argv[0]);
        return 1;
    }
    string txtfile(argv[1]);
    string vtkfile(argv[2]);
    int n = 0;
    if (argc == 4) {
        n = atoi(argv[3]);
    }

    Mesh mesh;
    mesh.convert_txt_2_vtk(txtfile, vtkfile, n);
    return 0;
}
