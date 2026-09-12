#include <CompactStar/Physics/Rotochemical/ScaledRKF45.hpp>
#include <iostream>
using namespace CompactStar::Physics::Rotochemical;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
int main(){try{
 require(SourceSHA256("")=="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855","SHA256 empty");
 require(SourceSHA256("abc")=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad","SHA256 abc");
 require(SourceSHA256(std::string(1000000,'a'))=="cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0","SHA256 multiblock");
 for(auto tol:{ComponentTolerances{},ComponentTolerances::Refined()}){
  std::unique_ptr<gsl_odeiv2_control,decltype(&gsl_odeiv2_control_free)> control(tol.Allocate(),gsl_odeiv2_control_free);
  for(size_t i=0;i<3;++i)for(double y:{0.,-1e-8,2.}){double level=0;require(gsl_odeiv2_control_errlevel(control.get(),y,987.,123.,i,&level)==GSL_SUCCESS,"GSL errlevel refusal");require(level==tol.absolute[i]+tol.relative*std::abs(y),"incorrect component tolerance mapping");}
 }
 for(double bad:{0.,-1.,double(INFINITY),double(NAN)}){ComponentTolerances tol;tol.relative=bad;bool failed=false;try{tol.Validate();}catch(...){failed=true;}require(failed,"invalid rtol accepted");tol={};tol.absolute[1]=bad;failed=false;try{tol.Validate();}catch(...){failed=true;}require(failed,"invalid component atol accepted");}
 bool null_spin=false;try{PrescribedDipoleHistory invalid(nullptr);}catch(...){null_spin=true;}require(null_spin,"missing owned spin token accepted");
 auto token=std::make_shared<RunDependencyToken>();PrescribedDipoleHistory history(token);token->generation++;bool stale_spin=false;try{history.Sample(0);}catch(...){stale_spin=true;}require(stale_spin,"stale spin callback accepted");
 std::cout<<"PASS SHA256 known answers and actual GSL component error levels; baseline/refined; invalid tolerances refused\n";return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
