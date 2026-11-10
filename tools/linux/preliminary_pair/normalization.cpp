// Export fixed geometry operands, not a new stellar background or trajectory.
#define NATIVE_STAGE_NO_MAIN
#include "../../../tests/native/stage.cpp"
#undef NATIVE_STAGE_NO_MAIN
#include <CompactStar/Physics/Driver/Thermal/Boundary/TbDefinition.hpp>
int main(int argc,char** argv){
 try{
  gsl_set_error_handler_off();require(argc==10,"profile thermal entry frozen coefficients archive sha work output");
  Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
  const auto ctx=pair.source->OrdinaryContext()->DriverContext();
  const auto index=CompactStar::Physics::Driver::Thermal::Boundary::FindTbIndex(*ctx.star,1e10);
  require(ctx.geo!=nullptr,"missing fixed geometry");
  std::ofstream f(argv[9]);f<<std::setprecision(17)<<"quantity\tvalue\n";
  f<<"R_km\t"<<ctx.geo->R()[-1]<<"\nM_km\t"<<ctx.geo->Mass()[-1]<<"\nexpnu_surface\t"<<ctx.geo->ExpNu()[-1]<<"\nexp2nu_surface\t"<<ctx.geo->Exp2Nu()[-1]<<"\nexpnu_boundary\t"<<ctx.geo->ExpNu()[index]<<"\nboundary_index\t"<<index<<'\n';
  f.flush();require(bool(f),"normalization export failed");return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
