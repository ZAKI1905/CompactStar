// ADR-0012 U9 independent flux and Stieltjes checks on exported A1 profiles.
#include "tests/rotation/hartle_monopole_reference.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
int main(int argc,char**argv){gsl_set_error_handler_off();std::cout<<std::setprecision(17);for(int arg=1;arg<argc;arg++){
 std::ifstream in(argv[arg]);std::string line;getline(in,line);hartle_ref::Background b;hartle_mono_ref::Background2 b2;std::vector<double> pm,pp,px;
 while(getline(in,line)){std::istringstream ss(line);double a[26];for(int i=0;i<26;i++)a[i]=0;for(int i=0;i<26&&ss;i++)ss>>a[i];
 b.r.push_back(a[0]);b.m.push_back(a[1]);b.p.push_back(a[3]);b.eps.push_back(a[2]);b.nu.push_back(a[5]);b.lambda.push_back(-.5*log(1-2*a[1]/a[0]));
 b2.r.push_back(a[0]);b2.m.push_back(a[1]);b2.p.push_back(a[3]);b2.eps.push_back(a[2]);b2.nu.push_back(a[5]);b2.nup.push_back(a[6]);b2.dedp.push_back(a[7]);b2.s.push_back(a[8]);b2.sp.push_back(a[9]);pm.push_back(a[10]);pp.push_back(a[11]);px.push_back(a[12]);}
 auto f=hartle_ref::Solve(b,1.,b.r.front());if(!f.ok)return 3;
 double R=b.R(),I=pow(R,4)*b2.sp.back()/6;std::cout<<argv[arg]<<" Iprod "<<I<<" Iflux "<<f.I_surface<<" Ivolume "<<f.I_volume;
 hartle_mono_ref::MHOptions opt;opt.I_exterior=I;hartle_mono_ref::StieltjesOptions st;st.refine=4;
 auto h=hartle_mono_ref::SolveStieltjes(b2,opt,st);if(!h.ok)return 4;
 std::cout<<" mh "<<h.mhat_R<<" ph "<<h.phat_R<<" xi "<<h.xihat_R<<" shell "<<h.shell_hat<<" deltaM "<<h.deltaM_hat;
 double err[3]={0,0,0};std::vector<double>*res[]={&h.mhat,&h.phat,&h.xihat};std::vector<double>*prod[]={&pm,&pp,&px};
 for(int k=0;k<3;k++){double peak=0;for(double x:*res[k])peak=std::max(peak,std::abs(x));for(size_t j=0;j<pm.size();j++)err[k]=std::max(err[k],std::abs((*res[k])[j]-(*prod[k])[j])/peak);std::cout<<" profile_norm"<<k<<" "<<err[k];}
 if(std::abs(I/f.I_surface-1)>1e-3)return 5;
 double prodDelta=pm.back()+4*M_PI*R*R*b2.eps.back()*px.back()+I*I/(R*R*R);
 std::cout<<" prod_deltaM "<<prodDelta<<" deltaM_rel "<<prodDelta/h.deltaM_hat-1;
 if(std::abs(prodDelta/h.deltaM_hat-1)>1e-4)return 6;
 for(double e:err)if(e>1e-4)return 7;
 std::cout<<" PASS\n";
}}
