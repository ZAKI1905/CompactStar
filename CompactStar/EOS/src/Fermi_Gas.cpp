/*
  Fermi_Gas Model Class
*/

#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/EOS/Fermi_Gas.hpp"
#include "CompactStar/EOS/Common.hpp"

using namespace Zaki::Physics ;
//==============================================================
//             Fermi_Gas Model Class
//==============================================================
// Constructor
CompactStar::Fermi_Gas::Fermi_Gas(const double& in_m) 
: Model({1e-12,10e+0}), m(in_m)
{
  SetName("Fermi_Gas") ;
}

//--------------------------------------------------------------
/// Destructor
CompactStar::Fermi_Gas::~Fermi_Gas() { }

//--------------------------------------------------------------
// Mass is in MeV
void CompactStar::Fermi_Gas::SetParams(const double& in_m)
{
  m = in_m ;
}

//--------------------------------------------------------------
// Energy density in g/cm^3 
// --> input rho has units of fm^{-3}
double CompactStar::Fermi_Gas::EDens(const double& in_rho) 
{

  // kF has units of fm^{-1}
  double kF = pow(3*M_PI*M_PI*in_rho, 1.0/3.0) ;

  // Converting the mass from MeV to fm^{-1}
  double m2 = m*MEV_2_INV_FM ;

  double mu = sqrt( m2*m2 + kF*kF ) ;

  // eps has unit of fm^{-4}
  double eps = mu*kF*(mu*mu - m2*m2/2) ;
        eps -= (pow(m2, 4)/2) * log((mu + kF) / m2) ;
        eps /= 4*M_PI*M_PI ;

  // converting eps unit from fm^{-4} to g/cm^3
  return eps*INV_FM4_2_G_CM3 ;
}

//--------------------------------------------------------------
// Pressure in units of dyne/cm^2 
//  for a given number density in_rho(1/fm^3)
double CompactStar::Fermi_Gas::Press(const double& in_rho) 
{
    // kF has units of fm^{-1}
  double kF = pow(3*M_PI*M_PI*in_rho, 1.0/3.0) ;

  // Converting the mass from MeV to fm^{-1}
  double m2 = m*MEV_2_INV_FM ;

  double mu = sqrt( m2*m2 + kF*kF ) ;

  // p has unit of fm^{-4}
  double p = mu*kF*(mu*mu - 5.*m2*m2/2) ;
        p += (3.*pow(m2, 4)/2) * log((mu + kF) / m2) ;
        p /= 12*M_PI*M_PI ;

  // converting p unit from fm^{-4} to g/cm^3
  return p*INV_FM4_2_Dyn_CM2 ;

}

//--------------------------------------------------------------

//==============================================================
