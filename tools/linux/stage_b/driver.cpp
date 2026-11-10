// Authorized Stage B scheduling wrapper; retained production methods and guards.
#define NATIVE_STAGE_NO_MAIN
#include "stage_b.inc"
#undef NATIVE_STAGE_NO_MAIN
int main(int argc,char** argv) {
 try {
  gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
  require(argc==13,"profile thermal entry frozen coefficients archive sha fresh-work output role reference-index tier-index");
  Inputs in;in.profile=argv[1];in.certificate=in.profile/"certificate.txt";in.thermal=argv[2];in.entry=argv[3];in.frozen=argv[4];in.coefficients=argv[5];in.work=argv[8];
  const std::string role=argv[10];const int reference=std::stoi(argv[11]),tier=std::stoi(argv[12]);
  require(role=="source"||role=="control","invalid role");
  require((reference>=0&&reference<5&&tier==-1)||(reference==-1&&tier>=0&&tier<3),"invalid independent task");
  auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);auto pair=NativePair(in,sequence);
  stage_role(role=="source"?pair.source:pair.control,2048*Year,argv[9],reference,tier);
  std::cout<<"Stage B task complete; accuracy assessment required\n";return 0;
 }catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
