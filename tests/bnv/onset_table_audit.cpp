// Validation only: original vs capped identities and strict serialized ordering.
#include "onset_capped_table.hpp"
#include <iostream>
int main(int argc,char**argv){try{
 structure1::require(argc==4,"out resolution cap");
 CompactStar::TrackRFreeGasThermodynamicProvider eos;
 auto out=std::filesystem::path(argv[1]);size_t n=std::stoull(argv[2]);int cap=std::stoi(argv[3]);
 auto c=onset_validation::generate(eos,out/"capped.tsv",n,1e-14,.616,cap);
 if(n<=8192){auto o=structure1::generate(eos,out/"original.tsv",n);
 std::ifstream a(c.file),b(o.file);std::string x((std::istreambuf_iterator<char>(a)),{}),y((std::istreambuf_iterator<char>(b)),{});
 structure1::require(x==y,"P5 identity failure");}
 std::cout<<"PASS rows "<<c.rows.size()<<" resolution "<<n<<" cap "<<cap<<'\n';
 }catch(const std::exception&e){std::cerr<<"REFUSED "<<e.what()<<'\n';return 1;}}
