// Linked production probe; all numbers emitted as exact binary64 inputs/outputs.
#include <CompactStar/EOS/TrackRFreeGasThermodynamics.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <fstream>
#include <string>
#include <set>
int main() {
  using namespace CompactStar;
  using namespace Zaki::Physics;
  TrackRFreeGasThermodynamicProvider p;
  auto e=ColdRelativisticIdealFermion::Electron();
  std::cout << std::hexfloat;
  std::cout << "constants\t" << NEUTRON_M_MEV << '\t' << PROTON_M_MEV << '\t' << ELECTRON_M_MEV << '\t' << MEV_2_INV_FM << '\n';
  std::cout << "onset\t" << p.NeutronOnsetBaryonDensityFm3() << '\n';
  std::cout << "muon_onset\t" << p.MuonOnsetBaryonDensityFm3() << '\n';
  const double mu=std::fma(NEUTRON_M_MEV-PROTON_M_MEV,NEUTRON_M_MEV+PROTON_M_MEV,ELECTRON_M_MEV*ELECTRON_M_MEV)/(2*NEUTRON_M_MEV);
  for(int offset : {-32,-8,-1,0,1,8,32}) {
    double x=mu;
    for(int i=0;i<std::abs(offset);++i) x=std::nextafter(x,offset<0?0.0:2.0);
    std::cout << "inverse\t" << x << '\t' << e.NumberDensityForChemicalPotentialFm3(x) << '\n';
  }
  for (double x : {0.8125,1.0625,1.75,2.125,8.5})
    std::cout << "inverse\t" << x << '\t' << e.NumberDensityForChemicalPotentialFm3(x) << '\n';
#ifdef __linux__
  std::ifstream maps("/proc/self/maps"); std::string line; std::set<std::string> paths;
  while(std::getline(maps,line)) { auto pos=line.find('/'); if(pos!=std::string::npos) paths.insert(line.substr(pos)); }
  for(const auto& path:paths) std::cout << "loaded\t" << path << '\n';
#endif
  for(auto fermion : {ColdRelativisticIdealFermion::Electron(),ColdRelativisticIdealFermion::Neutron()})
    for(double n : {1e-18,1e-12,1e-6,0.03125,0.5}) {
      auto v=fermion.Evaluate(n);
      std::cout << "local\t" << n << '\t' << fermion.RestMassEnergyMeV() << '\t' << v.fermi_momentum_MeV << '\t' << v.chemical_potential_MeV << '\t' << v.energy_density_MeV_fm3 << '\t' << *v.dchemical_potential_dn_MeV_fm3 << '\n';
    }
  std::cout << "types\t" << sizeof(float) << '\t' << std::numeric_limits<float>::digits << '\t' << sizeof(double) << '\t' << std::numeric_limits<double>::digits << '\t' << sizeof(long double) << '\t' << std::numeric_limits<long double>::digits << '\n';
}
