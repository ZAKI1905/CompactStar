// Read-only export of the retained thermal cache and zero-state normalizations.
#define NATIVE_STAGE_NO_MAIN
#include "../../../tests/native/stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
#include <CompactStar/Physics/BNV/CstarSegments.hpp>
int main(int argc,char** argv){
 try{
  gsl_set_error_handler_off();require(argc==10,"profile thermal entry frozen coefficients archive sha work output");
  Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
  std::filesystem::path out=argv[9];require(!std::filesystem::exists(out),"fresh export required");std::filesystem::create_directories(out);
  const BNV::CstarSegments thermal(pair.source->OrdinaryContext()->DriverContext());
  std::ofstream cells(out/"capacity.tsv");cells<<std::setprecision(17)<<"T_MeV\tC_erg_K\n";
  for(std::size_t i=0;i<thermal.Cells().temperature_MeV.size();++i)cells<<thermal.Cells().temperature_MeV[i]<<'\t'<<thermal.Cells().capacity_erg_K[i]<<'\n';
  ActualCheckpointContext context(pair.source,0);auto d=context.EvaluateDiagnostics(0,{{0,0,0}});
  std::ofstream f(out/"anchors.tsv");f<<std::setprecision(17)<<"quantity\tvalue\n";
  for(const auto& kv:DiagnosticValues(d))f<<kv.first<<'\t'<<kv.second<<'\n';
  f<<"sigma_e\t"<<d.sigma_count_s[0]<<"\nsigma_mu\t"<<d.sigma_count_s[1]<<'\n';
  cells.flush();f.flush();require(bool(cells)&&bool(f),"export failed");
  std::cout<<"PASS read-only thermal cache export; no integration\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
