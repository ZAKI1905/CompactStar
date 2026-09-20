#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
std::vector<std::string> Split(const std::string& line)
{std::istringstream in(line);std::vector<std::string> out;std::string value;while(std::getline(in,value,'\t'))out.push_back(value);return out;}
struct Table{std::vector<std::string> header;std::vector<std::vector<std::string>> rows;};
Table Read(const std::string& path)
{std::ifstream in(path);if(!in)throw std::runtime_error("missing trajectory table");Table out;std::string line;if(!std::getline(in,line))throw std::runtime_error("empty trajectory table");out.header=Split(line);while(std::getline(in,line))out.rows.push_back(Split(line));for(const auto& row:out.rows)if(row.size()!=out.header.size())throw std::runtime_error("trajectory table shape mismatch");return out;}
std::size_t Column(const Table& table,const std::string& name)
{for(std::size_t i=0;i<table.header.size();++i)if(table.header[i]==name)return i;throw std::runtime_error("missing trajectory column: "+name);}
double Value(const std::string& text){const double x=std::stod(text);if(!std::isfinite(x))throw std::runtime_error("nonfinite trajectory value");return x;}
}
int main(int argc,char** argv){try{
 if(argc!=3)throw std::runtime_error("usage: baseline.tsv refined.tsv");const auto base=Read(argv[1]),refined=Read(argv[2]);
 if(base.header!=refined.header||base.rows.size()!=refined.rows.size())throw std::runtime_error("ODE refinement schema/grid mismatch");
 const auto time=Column(base,"t_s");const std::array<std::string,3> names{{"x_state","eta_e_MeV","eta_mu_MeV"}};
 const std::array<double,3> atol{{1e-12,1e-18,1e-18}};std::array<double,3> maximum{};
 for(std::size_t row=0;row<base.rows.size();++row){if(base.rows[row][time]!=refined.rows[row][time])throw std::runtime_error("ODE refinement time grid mismatch");for(std::size_t i=0;i<3;++i){const auto c=Column(base,names[i]);const double a=Value(base.rows[row][c]),b=Value(refined.rows[row][c]);const double scaled=std::abs(a-b)/(atol[i]+1e-7*std::max(std::abs(a),std::abs(b)));maximum[i]=std::max(maximum[i],scaled);if(scaled>1)throw std::runtime_error("BA12 component-scaled difference failed: "+names[i]);}}
 std::cout.precision(17);std::cout<<"BA12_ODE_REFINEMENT PASS rows "<<base.rows.size()<<" max_x "<<maximum[0]<<" max_eta_e "<<maximum[1]<<" max_eta_mu "<<maximum[2]<<'\n';return 0;
 }catch(const std::exception& error){std::cerr<<"STOP "<<error.what()<<'\n';return 1;}}
