// Diagnostic-only frozen-Jacobian perturbation experiment; NOT the BNV RHS.
// No thermal/chemical energy ledger can be inferred from this surrogate.
#include <CompactStar/Physics/BNV/BreakpointStepper.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
using B=CompactStar::Physics::BNV::BreakpointStepper;
int main(int argc,char** argv) {
 try {
  if(argc!=2)throw std::runtime_error("expected matrix/state/cap input");
  gsl_set_error_handler_off();
  std::ifstream f(argv[1]);f.exceptions(std::ios::failbit|std::ios::badbit);
  double J[3][3]; B::State eq{};double cap;
  for(auto& r:J)for(auto& x:r)f>>x;for(auto& x:eq)f>>x;f>>cap;
  std::cout<<std::setprecision(17)<<"tier\taxis\tsign\tpolicy\tstep\ttime\tx\te\tm\trejected\trhs\twall_s\th_s\n";
  for(int tier=0;tier<3;++tier)for(int axis=0;axis<3;++axis)for(int sign:{-1,1})for(int policy=0;policy<4;++policy){
   double r=std::pow(10.,-7-2*tier);B::State a{{1e-12*std::pow(.01,tier),1e-18*std::pow(.01,tier),1e-18*std::pow(.01,tier)}};
   B::State y=eq;y[axis]+=sign*.25*(a[axis]+r*std::abs(eq[axis]));
   auto rhs=[&](double,const B::State& in,B::State& out,int){for(int i=0;i<3;++i){out[i]=0;for(int j=0;j<3;++j)out[i]+=J[i][j]*(in[j]-eq[j]);}};
   B solver(rhs,{},gsl_odeiv2_step_rkf45,r,a,policy==0?INFINITY:cap/std::pow(2.,policy-1));
   double t=0,h=1.0512181501887488e11;std::size_t count=0;
   auto start=std::chrono::steady_clock::now();
   const auto emit=[&](double dt){const auto s=solver.Stats();
    std::cout<<tier<<'\t'<<axis<<'\t'<<sign<<'\t'<<policy<<'\t'<<count<<'\t'<<t;
    for(double x:y)std::cout<<'\t'<<x;
    std::cout<<'\t'<<s.rejected<<'\t'<<s.rhs<<'\t'<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\t'<<dt<<'\n';};
   emit(0);
   while(t<5e11){if(++count>100000)throw std::runtime_error("work exhausted");auto s=solver.Advance(t,y,h,5e11);emit(s.right-s.left);}
  }
  std::cout.flush();if(!std::cout)throw std::runtime_error("output failure");
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
