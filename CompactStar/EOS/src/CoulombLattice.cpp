/*
  Coulomb Lattice Model Class
*/

// #include <gsl/gsl_math.h>

#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/EOS/CoulombLattice.hpp"
#include "CompactStar/EOS/Common.hpp"

using namespace Zaki::Physics ;
//==============================================================
//             Coulomb Lattice Model Class
//==============================================================
// Constructor
CompactStar::CoulombLattice::CoulombLattice() 
: Model({4.2144e-12,3.61317e-6}), 
  element("C", 6, 12)
{
  SetName("CoulombLattice") ;
}

//--------------------------------------------------------------
/// Destructor
CompactStar::CoulombLattice::~CoulombLattice() { }

//--------------------------------------------------------------
void CompactStar::CoulombLattice::SetElement(const Element& in_elem)
{
  element = in_elem ;
}

//--------------------------------------------------------------
// Energy density
double CompactStar::CoulombLattice::EDens(const double& in_rho) 
{
  double R = pow(3.*element.A/(4.*M_PI*in_rho), 1./3.);

  // element.mu() is the atomic mass in GeV
  double out = (element.mu() - element.Z*ELECTRON_M_GEV)*1e+3 
              - 0.9*element.Z*element.Z*Q_E_SQRD_MeVfm / R ;

  out *= in_rho/element.A ;

  // Converting the units
  out *= MEV_FM3_2_G_CM3 ;

  out += FermiEDens(pow(3.*M_PI*M_PI*in_rho*element.Z/element.A, 1./3.))*
          INV_FM4_2_G_CM3 ;

  return out ;
}
//--------------------------------------------------------------
// Pressure
double CompactStar::CoulombLattice::Press(const double& in_rho) 
{
  double out = -0.3*pow(4.*M_PI/3., 1./3.)
                *element.Z*element.Z*Q_E_SQRD_MeVfm
                *pow(in_rho/element.A, 4./3.) ;

  // Converting the units
  out *= MEV_FM3_2_Dyn_CM2 ;

  out += FermiPress(pow(3.*M_PI*M_PI*in_rho*element.Z/element.A, 1./3.))
        * INV_FM4_2_Dyn_CM2 ;

  return out ;
}

//--------------------------------------------------------------

//==============================================================
