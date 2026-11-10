// Selection/serialization tests only: no RHS, GSL, stellar or thermal solve.
#include "restart_schedule.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
struct Step {double t_left_s,t_right_s;std::array<double,3> y_right{};};
using History=std::vector<Step>;
void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
std::size_t select(const History& h,double a,double b) {
 return NativeQualification::RestartAtOrBeforeMidpoint(h,a,b)-h.cbegin();
}
void refuses(const History& h,double a,double b) {
 bool failed=false;try {select(h,a,b);}catch(const std::runtime_error&) {failed=true;}
 require(failed,"invalid coverage/history accepted");
}
int main(int argc,char** argv) {
 try {
  require(select({{0,1},{1,3},{3,10}},0,10)==1,"straddling terminal step");
  require(select({{0,1},{1,5},{5,10}},0,10)==1,"exact midpoint excluded");
  require(select({{0,1},{1,4},{4,6},{6,10}},0,10)==1,"chose endpoint after midpoint");
  require(select({{100,101},{101,104},{104,110}},100,110)==1,"absolute time origin");
  refuses({},0,10);refuses({{0,10}},0,10);refuses({{0,0},{0,10}},0,10);
  refuses({{0,1},{2,10}},0,10);refuses({{0,1},{1,9}},0,10);
  refuses({{0,1},{1,11}},0,10);refuses({{0,1},{1,10}},10,0);
  refuses({{0,1}},0,std::numeric_limits<double>::infinity());
  require(argc==1||argc==2,"optional accepted history TSV only");
  if(argc==2) {
   std::ifstream input(argv[1]);require(bool(input),"missing accepted history");
   std::string line;std::getline(input,line);History h;std::size_t count=0;
   while(std::getline(input,line)) {
    std::istringstream row(line);std::size_t ordinal,rejected;double step,next;Step s;
    row>>ordinal>>s.t_left_s>>s.t_right_s>>step>>next>>rejected;
    for(double& v:s.y_right)row>>v;
    require(bool(row)&&ordinal==++count&&rejected==0,"malformed retained history");
    h.push_back(s);
   }
   require(h.size()==12,"unexpected retained step count");
   const auto index=select(h,0,31557600.);
   require(index==10&&h[index].t_right_s==10025181.919443566,"wrong retained endpoint");
   require(31557600.-h[index].t_right_s>=31557600./2,"insufficient tail duration");
   std::ostringstream serialized;serialized<<std::setprecision(17)<<h[index].t_right_s;
   for(double v:h[index].y_right)serialized<<'\t'<<v;
   std::istringstream restored(serialized.str());double t;std::array<double,3> y;
   restored>>t;for(double& v:y)restored>>v;
   require(bool(restored)&&t==h[index].t_right_s&&y==h[index].y_right,"state roundtrip changed doubles");
   std::cout<<std::setprecision(17)<<"PASS retained ordinal="<<index+1<<" restart_s="<<t
            <<" tail_s="<<31557600.-t<<" exact state roundtrip\n";
  }
  std::cout<<"PASS 12 synthetic selection/refusal cases; no integration\n";
  return 0;
 } catch(const std::exception& e) {std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
