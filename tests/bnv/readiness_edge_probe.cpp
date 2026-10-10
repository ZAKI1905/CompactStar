// Preserve the accepted algorithm; characterize starts and safe refusals.
#include <CompactStar/Physics/BNV/BreakpointStepper.hpp>
#include <iostream>
#include <iomanip>
#include <string>
using B=CompactStar::Physics::BNV::BreakpointStepper;
void need(bool v,const char* s){if(!v)throw std::runtime_error(s);}
int main(){try{
 gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<"case\tmethod\tevents\tfinal_t\tresidual\tstatus\n";
 for(auto method:{gsl_odeiv2_step_rkf45,gsl_odeiv2_step_rk8pd})for(int direction:{-1,1})for(int side:{-1,0,1}){
  B solver([=](double,const B::State&,B::State& d,int){d={{direction*1e-12,2e-15,-3e-15}};},
    {-2,0,2},method,1e-12,{{1e-14,1e-18,1e-18}});
  double t=1e12,h=1;const double start=t,end=std::nextafter(std::nextafter(t,INFINITY),INFINITY);
  B::State y{{side*1e-17,0,0}};const auto initial=y;int n=0;
  while(t<end){need(++n<10,"micro-event loop");solver.Advance(t,y,h,end);}
  need(t==end,"endpoint lost");need(solver.Events().size()==std::size_t(side*direction<0),"unexpected event count");
  double residual=0;for(auto& e:solver.Events()){need(e.time_high>=e.time_low&&e.time>start,"invalid event bracket");residual=e.residual;}
  const double error=std::abs(y[0]-(initial[0]+direction*1e-12*(end-start)));need(error<1e-30,"state continuity lost");
  std::cout<<"start_"<<direction<<'_'<<side<<'\t'<<method->name<<'\t'<<solver.Events().size()<<'\t'<<t<<'\t'<<residual<<"\tPASS\n";
 }
 for(bool stationary:{false,true}){
  bool refused=false;std::string message;
  try{
   B solver([=](double t,const B::State&,B::State& d,int){d={{stationary?0.:7*std::cos(7*t),0,0}};},
      {-2,0,2},gsl_odeiv2_step_rkf45,10,{{10,10,10}});
   double t=0,h=1;B::State y{{stationary?0.:-.1,0,0}};solver.Advance(t,y,h,1);
  }catch(const std::exception& e){message=e.what();refused=message==(stationary?"unresolved stationary knot start":"unresolved reversing breakpoint trial");}
  need(refused,"expected specific fail-closed refusal missing");
  std::cout<<(stationary?"stationary":"turning")<<"\trkf45\t0\t0\t0\tPASS: "<<message<<'\n';
 }
 std::cout.flush();need(bool(std::cout),"output failure");
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
