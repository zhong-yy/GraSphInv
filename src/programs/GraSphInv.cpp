#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

#include "GaussNewtonInversion.h"
#include "timer.h"
class GraSphInv {
   public:
    GraSphInv();
    ~GraSphInv();

    void read_data_parameters(string data_para);
    void read_model_parameters(string model_para);
    void read_inversion_parameters(string inversion_para);
    void set_config_file(string file_name) { this->config_file = file_name; }
    void start_inversion();
    void write_result();

    int read_data_from_file(string file_name);
    istream& next_valid_line(istream& is, string& str);

    //处理注释，空格，和空行的函数
    // line，表示一行文本内容
    // comment_str，表示注释前导字符串，默认设置为#，也可以用//或者%
    // string删除字串，http://www.cplusplus.com/reference/string/string/erase/
    // string查找字串，http://www.cplusplus.com/reference/string/string/find_first_not_of/
    void line_process(std::string& line, const std::string comment_str = "#");

   protected:
    double height;
    double reference_level;
    int n_obs;
    Observation ob;
    VectorXd dobs;
    Mesh inv_mesh;
    GaussNewtonInversion* inv;

    int n_fields;
    unsigned long long field_flag;

    string data_file;
    string config_file;

    double min_size_r, min_size_lat, min_size_lon;
    double lat_model[2];
    double lon_model[2];
    double depth_model[2];
    double n_r, n_lat, n_lon;

    vector<double> noise_percentage;
    vector<double> equipment_noise;

    int GLQ_order;
    int integral_kernel;
    double as, ar, atheta, aphi, acrg;
    double depth_weighting_exponent;
    double start_lambda;
    double decreasing_rate;
    double target_misfit;  //数据误差收敛准则
    double cg_tol;         //共轭梯度收敛准则
    double stagnate_tol;
    int gauss_newton_iterations;
    int n_lambda;
    double cg_iteration_factor;
    double min_value, max_value;

    double refinement_tol;
    int interval_between_two_refinements;
    int max_refinement;
    int show_process_or_not;

    double relative_threshold;
    int method_id;

    string crg_model_file;
    string ref_model_file;

    int n_lat_crg_model, n_lon_crg_model, n_r_crg_model;
    int n_lat_pet_model, n_lon_pet_model, n_r_pet_model;

    bool use_crg;
    bool use_petr;

    string format_of_coordinates_crg;
    string format_of_coordinates_petr;
    int fast_dimension_crg;
    int fast_dimension_petr;

    string output_model_name;
};

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Usage: %s [configuration_file]\n", argv[0]);
        return 1;
    }

    assert(argc == 2);
    ifstream config(argv[1]);
    if (!config.good()) {
        cout << "Cannot read " << argv[1] << ", please check whether "
             << argv[1] << " exists" << endl;
        return 1;
    }

    cout << "Using inversion parameters from configuration file: " << argv[1]
         << endl
         << endl;

    GraSphInv inv_case;

    string data_para;
    string model_para;
    string inversion_para;

    string line;

    inv_case.next_valid_line(config, line);
    istringstream iss1(line);
    iss1 >> data_para;

    inv_case.next_valid_line(config, line);
    istringstream iss2(line);
    iss2 >> model_para;

    inv_case.next_valid_line(config, line);
    istringstream iss3(line);
    iss3 >> inversion_para;

    // this must called before reading model parameters
    inv_case.read_data_parameters(data_para);
    inv_case.read_model_parameters(model_para);
    inv_case.read_inversion_parameters(inversion_para);

    inv_case.start_inversion();
    inv_case.write_result();

    return 0;
}

GraSphInv::GraSphInv() {
    inv = NULL;
    reference_level = 6378137;
    use_crg = false;
    use_petr = false;
    format_of_coordinates_crg = "yxz";
    format_of_coordinates_petr = "yxz";
    fast_dimension_crg = 0;
    fast_dimension_petr = 0;
    as = 1;
    ar = 1;
    atheta = 1;
    aphi = 1;
    acrg = 1;
    method_id = 1;
    relative_threshold = 1;
}
GraSphInv::~GraSphInv() {
    if (inv != NULL) delete inv;
}

void GraSphInv::line_process(std::string& line, const std::string comment_str) {
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

istream& GraSphInv::next_valid_line(istream& input_stream, string& line) {
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

void GraSphInv::read_data_parameters(string data_para) {
    ifstream input_stream(data_para.c_str());

    string line;

    next_valid_line(input_stream, line);
    istringstream iss(line);
    iss >> data_file;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> height;
    cout << "Height of data from the reference surface is " << height << endl;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> reference_level;
    cout << "Reference level is " << setprecision(7) << reference_level << endl;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);

    iss >> n_fields;
    cout << endl;
    // cout << n_fields << " gravity component"<< ((n_fields == 1) ? (" is") :
    // ("s are")) << " used" << endl; dobs.resize(n_obs * n_fileds);

    field_flag = 0;
    vector<unsigned int> field_label(n_fields);
    next_valid_line(input_stream, line);
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

    // double noise_percentage, equipment_noise;
    // config >> noise_percentage >> equipment_noise;

    noise_percentage.resize(n_fields);

    equipment_noise.resize(n_fields);
    for (int i = 0; i < n_fields; i++) {
        next_valid_line(input_stream, line);
        iss.clear();
        iss.str("");
        iss.str(line);
        iss >> noise_percentage[i] >> equipment_noise[i];
    }

    read_data_from_file(data_file);
}

int GraSphInv::read_data_from_file(string data_file) {
    ifstream input_stream(data_file.c_str());
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

            r0 = reference_level + height;
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
    return n_obs;
}

void GraSphInv::read_model_parameters(string model_para) {
    ifstream input_stream(model_para.c_str());
    string line;
    next_valid_line(input_stream, line);
    istringstream iss1(line);
    iss1 >> lat_model[0] >> lat_model[1] >> n_lat;

    next_valid_line(input_stream, line);
    istringstream iss2(line);
    iss2 >> lon_model[0] >> lon_model[1] >> n_lon;

    next_valid_line(input_stream, line);
    istringstream iss3(line);
    iss3 >> depth_model[0] >> depth_model[1] >> n_r;
    cout << "Range of the latitude: ";
    cout << lat_model[0] << ", " << lat_model[1] << endl;
    cout << "Range of the longitude: ";
    cout << lon_model[0] << ", " << lon_model[1] << endl;
    cout << "Range of the depth: ";
    cout << depth_model[0] << ", " << depth_model[1] << endl;

    inv_mesh.generate_regular_geographic_mesh(
        lat_model, n_lat, lon_model, n_lon, depth_model, n_r, reference_level);

    inv_mesh.out_model_vtk("initial_mesh.vtk");
}

void GraSphInv::read_inversion_parameters(string inversion_para) {
    ifstream input_stream(inversion_para.c_str());

    string line;
    next_valid_line(input_stream, line);
    istringstream iss(line);
    iss >> GLQ_order;
    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> integral_kernel;  // 0: surface integral; otherwise: volume integral

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> method_id;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> relative_threshold;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> as >> ar >> atheta >> aphi >> acrg;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> depth_weighting_exponent;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> start_lambda >> n_lambda >> decreasing_rate;

    //   double final_lambda, final_misfit;
    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> target_misfit;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> stagnate_tol;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> cg_tol >> cg_iteration_factor;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> gauss_newton_iterations;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> min_value >> max_value;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> max_refinement;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> min_size_lat >> min_size_lon >> min_size_r;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> interval_between_two_refinements;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> refinement_tol;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> show_process_or_not;

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    int flag = 0;
    iss >> flag;
    if (flag == 0) {
        use_crg = false;
    } else {
        use_crg = true;
    }

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    if (use_crg) {
        iss >> crg_model_file >> format_of_coordinates_crg;
        if (format_of_coordinates_crg == "xyz") {
            iss >> n_lat_crg_model >> n_lon_crg_model >> n_r_crg_model;
        } else if (format_of_coordinates_crg == "yxz") {
            iss >> n_lon_crg_model >> n_lat_crg_model >> n_r_crg_model;
        } else if (format_of_coordinates_crg == "zxy") {
            iss >> n_r_crg_model >> n_lat_crg_model >> n_lon_crg_model;
        } else if (format_of_coordinates_crg == "zyx") {
            iss >> n_r_crg_model >> n_lon_crg_model >> n_lat_crg_model;
        } else if (format_of_coordinates_crg == "xzy") {
            iss >> n_lat_crg_model >> n_r_crg_model >> n_lon_crg_model;
        } else if (format_of_coordinates_crg == "yzx") {
            iss >> n_lon_crg_model >> n_r_crg_model >> n_lat_crg_model;
        }
        iss >> fast_dimension_crg;
    }

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    flag = 0;
    iss >> flag;
    if (flag == 0) {
        use_petr = false;
    } else {
        use_petr = true;
    }

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    if (use_petr) {
        iss >> ref_model_file >> format_of_coordinates_petr;
        if (format_of_coordinates_petr == "xyz") {
            iss >> n_lat_pet_model >> n_lon_pet_model >> n_r_pet_model;
        } else if (format_of_coordinates_petr == "yxz") {
            iss >> n_lon_pet_model >> n_lat_pet_model >> n_r_pet_model;
        } else if (format_of_coordinates_petr == "zxy") {
            iss >> n_r_pet_model >> n_lat_pet_model >> n_lon_pet_model;
        } else if (format_of_coordinates_petr == "zyx") {
            iss >> n_r_pet_model >> n_lon_pet_model >> n_lat_pet_model;
        } else if (format_of_coordinates_petr == "xzy") {
            iss >> n_lat_pet_model >> n_r_pet_model >> n_lon_pet_model;
        } else if (format_of_coordinates_petr == "yzx") {
            iss >> n_lon_pet_model >> n_r_pet_model >> n_lat_pet_model;
        }
        iss >> fast_dimension_petr;
    }

    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> output_model_name;
}
void GraSphInv::start_inversion() {
    string line;

    Timer timer;
    timer.start();
    inv = new GaussNewtonInversion(inv_mesh, ob, field_flag);
    inv->set_GLQ_order(GLQ_order);
    inv->set_dobs(dobs, noise_percentage, equipment_noise);
    inv->set_depth_weighting(depth_weighting_exponent);
    inv->set_weights_of_objectives(as, ar, atheta, aphi, acrg);
    inv->set_max_lambda(start_lambda);
    inv->set_max_GN_iterations(gauss_newton_iterations);
    inv->set_lambda_decreasing_rate(decreasing_rate);
    inv->set_n_lambda(n_lambda);
    inv->set_CG_parameter(cg_tol, cg_iteration_factor);
    inv->set_interval_between_refinements(interval_between_two_refinements);
    inv->set_refinement_percentage(refinement_tol);
    inv->set_max_refinement_number(max_refinement);
    inv->set_min_cell_size_in_adaptive_mesh(min_size_lat, min_size_lon,
                                            min_size_r);
    inv->set_stagnation_tolerance(stagnate_tol);
    inv->set_target_misfit(target_misfit);

    if (use_crg) {
        inv->create_crg_model_from_data(
            crg_model_file, n_lat_crg_model, n_lon_crg_model, n_r_crg_model,
            format_of_coordinates_crg, fast_dimension_crg);
    }
    if (use_petr) {
        inv->create_ref_model_from_data(
            ref_model_file, n_lat_pet_model, n_lon_pet_model, n_r_pet_model,
            format_of_coordinates_petr, fast_dimension_petr);
    }

    // if()

    inv->set_integral_kernel_type(integral_kernel);
    inv->set_method_id(this->method_id);
    inv->set_compression_threshold(relative_threshold);
    // inv->set_use_wavelet(true);
    // inv->set_compression_threshold(0.01);

    int Nm = inv_mesh.n_elems();
    VectorXd m_min = VectorXd::Constant(Nm, min_value);
    VectorXd m_max = VectorXd::Constant(Nm, max_value);
    VectorXd m_ini = VectorXd::Constant(Nm, 0);

    inv->set_min_max(m_min, m_max);
    inv->set_m_ini(m_ini);

    cout << endl;

    cout << "Start inversion" << endl;

    double computation_time;

    VectorXd m_result, d_pre;

    if (show_process_or_not != 0) {
        // cout << "The inversion process will be recorded" << endl;
        inv->record_every_iteration();
    } else {
        // cout << "The inversion process will not be recorded" << endl;
    }

    inv->display_inversion_parameters();

    inv->invert();
    timer.stop();
    computation_time = timer.getElapsedTimeInSec();
    ofstream out_info("info");
    out_info << "1. Initial inverion mesh" << endl;
    out_info << "Range of the latitude: ";
    out_info << lat_model[0] << ", " << lat_model[1] << endl;
    out_info << "Range of the longitude: ";
    out_info << lon_model[0] << ", " << lon_model[1] << endl;
    out_info << "Range of the depth: ";
    out_info << depth_model[0] << ", " << depth_model[1] << endl;
    out_info << "Number of elements along radius in the initial mesh " << n_r
             << endl;
    out_info << "Number of elements along latitude in the initial mesh "
             << n_lat << endl;
    out_info << "Number of elements along longitude in the initial mesh "
             << n_lon << endl;
    out_info << "Total number of elements in the initial mesh "
             << n_r * n_lat * n_lon << endl
             << endl;
    out_info << "2. Final mesh" << endl;
    out_info << "Element number " << inv->number_unknowns() << endl << endl;

    out_info << "3. Evaluation of the inversion" << endl;
    out_info << "Final misfit " << inv->get_final_misfit() << endl;
    out_info << "Running time " << computation_time << " seconds" << endl;

    cout << "Running time of inversion:" << computation_time << "s" << endl;

    cout << "Inversion completed" << endl;
}

void GraSphInv::write_result() {
    cout << "Writing the inversion model into file ..." << endl;
    inv->result2vtk(output_model_name);
    inv->result2txt(output_model_name);
    inv->result2netcdf(output_model_name);

    cout << "Writing predicted data ..." << endl;
    inv->output_obs_data("dobs");
    inv->output_predicted_data("dpredicted");
}
