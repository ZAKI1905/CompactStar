#include <algorithm>
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
constexpr double Year=31557600.0,KbErg=1.380649e-16,KbMeV=8.617333262145e-11;
std::vector<std::string> Split(const std::string& line){std::istringstream in(line);std::vector<std::string> out;std::string x;while(std::getline(in,x,'\t'))out.push_back(x);return out;}
struct Table{std::vector<std::string> h;std::vector<std::vector<std::string>> r;};
Table Read(const std::string& path){std::ifstream in(path);if(!in)throw std::runtime_error("missing table");Table t;std::string line;if(!std::getline(in,line))throw std::runtime_error("empty table");t.h=Split(line);while(std::getline(in,line))t.r.push_back(Split(line));for(const auto& r:t.r)if(r.size()!=t.h.size())throw std::runtime_error("table shape mismatch");return t;}
std::size_t C(const Table&t,const std::string& n){for(std::size_t i=0;i<t.h.size();++i)if(t.h[i]==n)return i;throw std::runtime_error("missing column: "+n);}
double V(const Table&t,const std::vector<std::string>&r,const std::string&n){const double x=std::stod(r[C(t,n)]);if(!std::isfinite(x))throw std::runtime_error("nonfinite value: "+n);return x;}
double HM(double x){const double u2=(x/M_PI)*(x/M_PI);return x/(M_PI*M_PI)*(14680+u2*(7560+u2*(840+24*u2)))/11513;}
double D_HM(double x){const double u2=(x/M_PI)*(x/M_PI);return (14680+3*7560*u2+5*840*u2*u2+7*24*u2*u2*u2)/(11513*M_PI*M_PI);}
double Increment(double x){const double u2=(x/M_PI)*(x/M_PI);return u2*(22020+u2*(5670+u2*(420+9*u2)))/11513;}
}
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("usage: qss.tsv coefficients.tsv");const auto q=Read(argv[1]),coeff=Read(argv[2]);if(coeff.r.empty())throw std::runtime_error("empty coefficient table");
 const auto& c=coeff.r.front();const std::array<std::array<double,2>,2> z{{{{V(coeff,c,"Z00"),V(coeff,c,"Z01")}},{{V(coeff,c,"Z10"),V(coeff,c,"Z11")}}}};
 const std::array<double,2> L{{V(coeff,c,"Ltilde_Me"),V(coeff,c,"Ltilde_Mmu")}};std::array<double,2> max_balance{},max_tau_ratio{};double worst_bound=-1e300;std::size_t samples=0;
 for(const auto& row:q.r){const double t=V(q,row,"t_s");const std::array<double,2> xi{{V(q,row,"xi_e"),V(q,row,"xi_mu")}},eta{{V(q,row,"eta_e_MeV"),V(q,row,"eta_mu_MeV")}},R{{V(q,row,"R_e_count_s"),V(q,row,"R_mu_count_s")}},sigma{{V(q,row,"sigma_e_count_s"),V(q,row,"sigma_mu_count_s")}};const double T=V(q,row,"Tinf_K");
   for(double x:xi){const double bound=Increment(x)-x*HM(x);worst_bound=std::max(worst_bound,bound);if(bound>.467659+1e-10)throw std::runtime_error("BA17 modified-Urca bound failed");}
   if(t<4e5*Year)continue;++samples;std::array<double,2>d{};for(std::size_t i=0;i<2;++i){d[i]=L[i]/KbErg*std::pow(T,7)*D_HM(xi[i])/(KbMeV*T);const double resolution=std::abs(d[i])*(1e-18+1e-7*std::abs(eta[i]));max_balance[i]=std::max(max_balance[i],std::abs(R[i]+sigma[i])/std::max(std::abs(sigma[i]),resolution));}
   const double a=z[0][0]*d[0],b=z[0][1]*d[1],cc=z[1][0]*d[0],dd=z[1][1]*d[1],trace=a+dd,det=a*dd-b*cc,disc=trace*trace-4*det;if(!(disc>=0))throw std::runtime_error("complex QSS modes");const double slow=.5*(trace-std::sqrt(disc));if(!(slow>0))throw std::runtime_error("nonpositive QSS mode");const double ratio=(1/slow)/t;max_tau_ratio[0]=std::max(max_tau_ratio[0],ratio);max_tau_ratio[1]=std::max(max_tau_ratio[1],ratio);
 }
 if(samples==0)throw std::runtime_error("empty QSS terminal interval");for(std::size_t i=0;i<2;++i)if(max_balance[i]>.05||max_tau_ratio[i]>.10)throw std::runtime_error("BA16 linear-QSS criterion failed");
 std::cout.precision(17);std::cout<<"BA16_QSS PASS samples "<<samples<<" balance_e "<<max_balance[0]<<" balance_mu "<<max_balance[1]<<" tau_over_elapsed "<<max_tau_ratio[0]<<"\nBA17_BOUND PASS worst "<<worst_bound<<'\n';return 0;
 }catch(const std::exception&e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
