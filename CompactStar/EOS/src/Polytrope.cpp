/*
  Polytrope Model Class
*/

#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/EOS/Polytrope.hpp"
#include "CompactStar/EOS/Common.hpp"

using namespace Zaki::Physics ;
//==============================================================
//             Polytrope Model Class
//==============================================================
// Constructor
CompactStar::Polytrope::Polytrope() 
: Model({1e-12,1e-6})
{
  SetName("Polytrope") ;
}

//--------------------------------------------------------------
/// Destructor
CompactStar::Polytrope::~Polytrope() { }

//--------------------------------------------------------------
void CompactStar::Polytrope::SetParams(const double& in_K, 
                                       const double& in_gamma)
{
  K = in_K ;
  gamma = in_gamma ;
}

//--------------------------------------------------------------
// Energy density
double CompactStar::Polytrope::EDens(const double& in_rho) 
{
  return in_rho*938.91897*MEV_FM3_2_G_CM3 ;
}

//--------------------------------------------------------------
// Pressure for a given number density in_rho(1/fm^3)
double CompactStar::Polytrope::Press(const double& in_rho) 
{
  return K*pow(in_rho*938.91897*MEV_2_INV_FM, gamma)*INV_FM4_2_Dyn_CM2 ;
}

//--------------------------------------------------------------
/// Inverse function to "EDens"
double CompactStar::Polytrope::GetRho(const double& in_e) const
{
  return in_e/938.91897/MEV_FM3_2_G_CM3 ;
}
//==============================================================
