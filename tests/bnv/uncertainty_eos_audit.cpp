// Read-only diagnostic of the unchanged generator's strict ordering refusal.
#include "../eos/structure1/table.hpp"
#include <iostream>
int main(){
 using namespace structure1;TrackRFreeGasThermodynamicProvider eos;
 std::set<double> grid;const size_t resolution=32768;const double lower=1e-14,upper=.616;
 for(size_t i=0;i<=resolution;++i)grid.insert(lower*std::pow(upper/lower,double(i)/resolution));
 grid.insert(lower);grid.insert(upper);
 for(double t:{eos.NeutronOnsetBaryonDensityFm3(),eos.MuonOnsetBaryonDensityFm3()}){
  size_t per=resolution/128;int decades=t==eos.NeutronOnsetBaryonDensityFm3()?8+int(std::log2(double(resolution)/1024)):8;
  grid.insert(t);for(size_t i=per;i<=size_t(decades)*per;++i){double d=std::pow(10.,-double(i)/per);for(double n:{t*(1-d),t*(1+d)})if(n>lower&&n<upper)grid.insert(n);}
 }
 std::cout<<std::setprecision(17)<<"n_previous\tn_current\tP_previous\tP_current\te_previous\te_current\n";
 FreeGasBarotropePoint last{};bool first=true;unsigned count=0;
 for(double n:grid){if(n>upper)continue;auto p=eos.BarotropeAt(n);
  if(!first&&!(p.pressure_MeV_fm3>last.pressure_MeV_fm3&&p.energy_density_MeV_fm3>last.energy_density_MeV_fm3)){
   std::cout<<last.n_B_fm3<<'\t'<<n<<'\t'<<last.pressure_MeV_fm3<<'\t'<<p.pressure_MeV_fm3<<'\t'<<last.energy_density_MeV_fm3<<'\t'<<p.energy_density_MeV_fm3<<'\n';++count;
  }last=p;first=false;
 }
 std::cerr<<"ordering_failures "<<count<<'\n';return count?0:1;
}
