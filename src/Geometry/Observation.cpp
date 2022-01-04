#include "Observation.h"

Observation::Observation() : obs(0) {
  n_obs = 0;
}

Observation::~Observation() {}

void Observation::read_site(const string& name) {
  ifstream is;
  is.open(name.c_str());
  assert(is.good());
  is >> n_obs;
  cout << "\nReading observation sites ...\nThe total number of sites is "
       << n_obs << ".\n";
  unsigned lab;
  obs.resize(n_obs);
  double x, y, z;
  for (unsigned i = 0; i < n_obs; i++) {
    assert(is.good());
    is >> lab;
    // assert(lab == i + 1);
    is >> x >> y >> z;
    obs[i].set_xyz(x, y, z);
  }
  cout << "Reading sites done!\n\n";
}

void Observation::read_geographic_points(const string& name) {
  ifstream is;
  is.open(name.c_str());
  assert(is.good());
  string line;

  while (std::getline(is, line)) {
    line_process(line, "#");
    if (line.empty())
      continue;
    std::istringstream iss(line);
    iss >> n_obs;
    break;
  }

  cout << "\nReading observation sites ...\nThe total number of sites is "
       << n_obs << ".\n";
  unsigned lab;
  obs.resize(n_obs);
  double longitude, latitude, r;
  for (unsigned i = 0; i < n_obs; i++) {
    assert(is.good());
    while (std::getline(is, line)) {
      line_process(line, "#");
      if (line.empty())
        continue;
      std::istringstream iss(line);
      iss >> lab;
      iss >> longitude >> latitude >>
          r;  // read longitude, latitude in degree and radius in meter
      break;
    }
    double phi, theta;
    theta = (90.0 - latitude) * GS::PI /
            180.0;  // transform latitude in degree to phi in radian
    phi = (longitude + 180.0) * GS::PI /
          180.0;  // transform longitude in degree to phi in radian

    obs[i].set_xyz(r, theta, phi);
  }
  cout << "Reading sites done!\n\n";
}

void Observation::add_point(const Point& p) {
  obs.push_back(p);
  n_obs++;
}
void Observation::add_point(double x, double y, double z) {
  obs.push_back(Point(x, y, z));
  n_obs++;
}
void Observation::add_point_pre(double x, double y, double z) {
  obs.insert(obs.begin(), Point(x, y, z));
  n_obs++;
}

void Observation::generate(int nx,
                           double start_x,
                           double end_x,
                           int ny,
                           double start_y,
                           double end_y,
                           double z0) {
  obs.clear();
  obs.resize(nx * ny);
  double xspace = (end_x - start_x) / (nx - 1);
  double yspace = (end_y - start_y) / (ny - 1);
  for (int i_y = 0; i_y < ny; i_y++) {
    for (int i_x = 0; i_x < nx; i_x++) {
      obs[i_y * nx + i_x].set_xyz(start_x + i_x * xspace,
                                  start_y + i_y * yspace, z0);
    }
  }
  this->n_obs = nx * ny;
}
void Observation::generate_geographic(int n_r,
                                      double start_r,
                                      double end_r,
                                      int n_lat,
                                      double start_lat,
                                      double end_lat,
                                      int n_lon,
                                      double start_lon,
                                      double end_lon) {
  double start_theta = (90.0 - start_lat) * GS::PI / 180.0;
  double end_theta = (90.0 - end_lat) * GS::PI / 180.0;
  double start_phi = (180 + start_lon) * GS::PI / 180.0;
  double end_phi = (180 + end_lon) * GS::PI / 180.0;

  this->generate(n_r, start_r, end_r, n_lat, start_theta, end_theta, n_lon,
                 start_phi, end_phi);
}
void Observation::generate(int n_r,
                           double start_r,
                           double end_r,
                           int n_theta,
                           double start_theta,
                           double end_theta,
                           int n_phi,
                           double start_phi,
                           double end_phi) {
  obs.clear();
  obs.resize(n_r * n_theta * n_phi);
  double r_space;
  if (n_r == 1) {
    r_space = 0;
  } else {
    r_space = (end_r - start_r) / (n_r - 1.0);
  }
  double theta_space = (end_theta - start_theta) / (n_theta - 1.0);
  double phi_space = (end_phi - start_phi) / (n_phi - 1.0);
  for (int i_r = 0; i_r < n_r; i_r++) {
    for (int i_phi = 0; i_phi < n_phi; i_phi++) {
      for (int i_theta = 0; i_theta < n_theta; i_theta++) {
        obs[i_r * (n_theta * n_phi) + i_phi*n_theta+i_theta].set_xyz(
            start_r + i_r * r_space, start_theta + i_theta * theta_space,
            start_phi + i_phi * phi_space);
          }
    }
  }
  this->n_obs = n_r * n_theta * n_phi;
}

ostream& operator<<(ostream& output, Observation& observation) {
  output << "#Number of observation points: " << observation.obs.size() << '\n';
  output << "#r\ttheta\tphi\n";
  for (int i = 0; i < observation.n_obs; i++) {
    output << '\n';
    double* xyz = observation.obs[i].get_xyz();
    output << left << setw(15) << xyz[0] << left << setw(15)
           << (90 - xyz[1] * 180.0 / GS::PI) << left << setw(15)
           << xyz[2] * 180.0 / GS::PI - 180.0;
  }
  return output;
}
