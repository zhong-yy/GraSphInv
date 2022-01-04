#ifndef _OBSERVATION
#define _OBSERVATION
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "gs.h"
#include "point.h"
using namespace std;
class Observation {
 public:
  Observation();
  ~Observation();
  void add_point(const Point& p);
  void add_point(double x, double y, double z);
  void add_point_pre(double x, double y, double z);

  void generate_geographic(int n_r,
                           double start_r,
                           double end_r,
                           int n_lat,
                           double start_lat,
                           double end_lat,
                           int n_lon,
                           double start_lon,
                           double end_lon);
  void generate(int n_r,
                double start_r,
                double end_r,
                int n_theta,
                double start_theta,
                double end_theta,
                int n_phi,
                double start_phi,
                double end_phi);

  void generate(int nx,
                double start_x,
                double end_x,
                int ny,
                double start_y,
                double end_y,
                double z0);

  void read_site(const string& name);
  void read_geographic_points(const string& name);
  const Point& operator()(unsigned i) const {
    assert(i < n_obs);
    return obs[i];
  }
  void set_obs(unsigned i, double x, double y, double z) {
    assert(i < n_obs);
    this->obs[i].set_xyz(x, y, z);
  }
  unsigned int get_n_obs() const { return n_obs; }
  friend ostream& operator<<(ostream& output, Observation& observation);

  void clear(){
    this->obs.clear();
    this->n_obs=0;
  }

 private:
  vector<Point> obs;
  unsigned int n_obs;

  void line_process(string& line, const string comment_str = "#") {
    // if(line.size()==1&&(line[0]=='\n'||line[0]=='\r')){
    // 	cout<<"23333"<<endl;
    // }
    for (char& c : line)
      if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n')
        c = ' ';
    line.erase(0, line.find_first_not_of(" "));
    line.erase(line.find_last_not_of(" ") + 1);
    int n_comment_start = line.find_first_of(comment_str);
    if (n_comment_start != string::npos)
      line.erase(n_comment_start);
  }
};
#endif