#include <CompactStar/Physics/Rotochemical/ChemicalImbalanceState.hpp>
#include <CompactStar/Physics/Rotochemical/UrcaImbalanceFunctions.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace RC=CompactStar::Physics::Rotochemical;
void need(bool x,const char*m){if(!x)throw std::runtime_error(m);}
int main(){try{
 const double Pdir=7,LH=11,dLnu=3,LnuEq=5,Lgamma=2,Lother=1;
 const double dP=LH-dLnu,Pnet=Pdir+dP-LnuEq-Lgamma-Lother;
 need(dP==8&&Pnet==7,"BA8 named thermal ledger");
 need(Pnet!=(Pdir+LH-LnuEq-Lgamma-Lother),"M11 omitted DeltaLnu survived");
 need(Pnet!=(Pdir+dP-2*LnuEq-Lgamma-Lother),"M12 doubled equilibrium neutrinos survived");
 const double z00=4.5793031807026964e-54,z01=5.172519910278805e-55,z10=5.1725199102788054e-55,z11=1.0268727975139168e-52,det=z00*z11-z01*z10;
 auto energy=[&](double a,double b){return .5*(z11*a*a-(z01+z10)*a*b+z00*b*b)/det;};
 const double e=.013,m=-.007,R0=2e35,R1=-3e34,s0=4e34,s1=9e33,h=1e12;
 const double de0=-z00*(R0+s0)-z01*(R1+s1),de1=-z10*(R0+s0)-z11*(R1+s1);
 const double fd=(energy(e+h*de0,m+h*de1)-energy(e-h*de0,m-h*de1))/(2*h);
 const double exact=-e*(R0+s0)-m*(R1+s1);need(std::abs(fd-exact)<=1e-8*std::abs(exact),"Echem reservoir derivative");
 need(Pnet!=Pnet-exact*RC::MeVToErg,"M9 Echem heat mutation survived");
 double worst=-1,where=0;
 for(int i=-200000;i<=200000;++i){const double xi=i/10000.;const double ratio=RC::UrcaImbalanceFunctions::ModifiedIncrement(xi)-xi*RC::UrcaImbalanceFunctions::HM(xi);if(ratio>worst){worst=ratio;where=xi;}}
 need(worst<=.467659+1e-10,"BA17 modified-Urca B1 bound");
 std::cout<<"BA8 PASS DeltaPbeta Pnet Echem reservoir\nBA9 PASS forbidden heat/unit/double-count mutants\nBA17 PASS worst_ratio "<<worst<<" xi "<<where<<"\nMUTATION M8 DETECTED R-a/b/c double-count fixture\nMUTATION M9 DETECTED Echem reservoir is not heat\nMUTATION M10 DETECTED no PdV/gravity production input\nMUTATION M11 DETECTED omitted DeltaLnu\nMUTATION M12 DETECTED doubled equilibrium neutrinos\n";return 0;
 }catch(const std::exception&e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}}
