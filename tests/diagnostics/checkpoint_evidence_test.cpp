#define main retained_checkpoint_contract_main
#include "../bnv/passive_checkpoint_output_contract.cpp"
#undef main
#include "checkpoint_evidence.hpp"
#include <sstream>
int main(int argc,char** argv){try{
 Require(argc==2,"fresh test output required");const std::filesystem::path out=argv[1];Require(!std::filesystem::exists(out),"test output exists");
 auto audit=std::make_shared<FactoryAudit>();audit->equilibrium_energy_test=true;audit->derivative={{1,0,0}};
 BNV::Rk8pdCheckpointReconstructor solver;auto q=solver.Reconstruct(Interior(),{Step()},Factory(audit));
 Require(q.status==BNV::CheckpointStatus::Qualified,"synthetic state unexpectedly refused");solver.EvaluateDiagnostics(q,Factory(audit));
 Require(!q.diagnostics_evaluated&&q.diagnostic_U_O.at("Echem")>.20*q.diagnostic_F_P.at("Echem"),"deliberate Echem diagnostic failure absent");
 auto before=q;CheckpointEvidence::Write(out/"diagnostic",q,Step(),BNV::Rkf45Configuration{});
 Require(q.status==before.status&&q.O2_state==before.O2_state&&q.diagnostic_U_O==before.diagnostic_U_O&&q.diagnostic_F_P==before.diagnostic_F_P,"report changed acceptance or values");
 std::ifstream f(out/"diagnostic/diagnostics.tsv");std::string line;bool found=false;
 while(std::getline(f,line)){if(line.rfind("Echem\t",0)!=0)continue;std::stringstream s(line);std::string field;std::vector<std::string> row;while(std::getline(s,field,'\t'))row.push_back(field);
 Require(row[2]=="MeV"&&std::stod(row[5])==q.diagnostic_U_O.at("Echem")&&std::stod(row[6])==q.diagnostic_F_P.at("Echem")&&std::stod(row[7])==.20*q.diagnostic_F_P.at("Echem")&&row[15]=="0","failure fields were not preserved exactly");found=true;}
 Require(found,"failed channel missing");
 audit->varying_derivative=true;audit->derivative={{1e-4,1e-4,1e-4}};auto state=solver.Reconstruct(Interior(),{Step()},Factory(audit));
 Require(state.status==BNV::CheckpointStatus::NumericallyUnresolved,"synthetic state failure absent");CheckpointEvidence::Write(out/"state",state,Step(),BNV::Rkf45Configuration{});
 std::cout<<"PASS deliberate diagnostic/state failures persisted with exact U/F/units/time; reporting preserves status and numerical values\n";
 return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
