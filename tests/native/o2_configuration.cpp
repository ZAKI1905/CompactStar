// Configuration contract for the owner's narrowly scoped Stage A recovery.
#include <CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp>
#include <iostream>
#include <iomanip>
#include <stdexcept>
namespace BNV=CompactStar::Physics::BNV;
int main() {
 try {
  BNV::Rk8pdCheckpointReconstructor r;
  if(r.O1().relative!=1e-12||r.O1().absolute!=BNV::CheckpointState{{1e-17,1e-23,1e-23}}
    ||r.O2().relative!=1e-13||r.O2().absolute!=BNV::CheckpointState{{1e-18,5e-27,5e-27}})
    throw std::runtime_error("O1/O2 configuration differs from sole approved change");
  std::cout<<std::setprecision(17)<<"PASS O1 unchanged; O2 relative/thermal unchanged; chemical absolute="
           <<r.O2().absolute[1]<<','<<r.O2().absolute[2]<<" MeV\n";
  return 0;
 } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
