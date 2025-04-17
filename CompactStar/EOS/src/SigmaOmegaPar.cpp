/*
  Sigma-Omega Param Struct
*/

#include "iostream"

#include <Zaki/Physics/Constants.hpp>
// #include <Zaki/Math/GSLFuncWrapper.hpp>

#include "CompactStar/EOS/SigmaOmegaPar.hpp"
// #include "CompactStar/EOS/Common.hpp"

using namespace Zaki::Physics ;
//==============================================================
//                   Sigma-Omega Par Struct
//==============================================================
CompactStar::SigmaOmegaPar& 
CompactStar::SigmaOmegaPar::Derive
(const CompactStar::NuclMatterPar& n_par)
{
  //.........................................................
  //                    Initiliazing 
  //.........................................................

  // Fermi momentum: k            in fm^-1
  double k = n_par.kf() ;

  // Effective mass: m*           in fm^-1
  double ms = n_par.eff_m*MN ; 

  // E(k) = sqrt [ k^2 + m*^2 ]   in fm^-1
  double e = sqrt(k*k + ms*ms) ;

  // Binding energy               in fm^-1
  double eb = n_par.bind_e*Zaki::Physics::MEV_2_INV_FM ;

  // Saturation density
  double rho = n_par.sat_den ;

  // Compression factor           in fm^-1
  double comp = n_par.compr*Zaki::Physics::MEV_2_INV_FM ;

  // Symmetry energy              in fm^-1
  double a_sym_fm = n_par.a_sym*Zaki::Physics::MEV_2_INV_FM ;
  //.........................................................

  // First parameter is found here:
  gome_m_sqrd = (MN + eb -  e) / rho ;

  // In case of rho meson:
  grho_m_sqrd = 8./rho*(a_sym_fm - k*k/e/6.) ;
  //..........................................................
  //                        Intergals
  //..........................................................
  double x  = k / ms ; double t = sqrt(1 + x*x) ;
  double I1 = 2./pow(M_PI,2) * pow(ms,2)
              *(  x*t/2 + x/t - 1.5 * log(x+t) ) ;
  double I2 = 0.5/pow(M_PI,2) * pow(ms,4)
              *(  x*pow(t,3) - x*t/2 - 0.5 * log(x+t) ) ;
  double I3 = 1./pow(M_PI,2) * pow(ms,3)
              *(  x*t - log(x+t) ) ;
  //..........................................................

  //..........................................................
  // Equation coefficients in Eq. (4.187) on P. 160
  //..........................................................
  double alp_1 = comp 
                - gome_m_sqrd * 6. * pow(k, 3)/pow(M_PI,2)
                - 3. * k * k / e ;
  double bet_1 = 2 * (MN - ms) * alp_1 ;
  double gam_1 = 3 * pow(MN - ms, 2) * alp_1 ;
  double del_1 = -alp_1*I1 
                - 6. * pow(k, 3) / pow(M_PI,2) * pow(ms/e,2) ;

  double alp_2 = 0.5 * pow(MN - ms, 2) ;
  double bet_2 = 1./3. * MN * pow(MN - ms, 3) ;
  double gam_2 = 0.25 * pow(MN - ms, 4) ;
  double del_2 = rho * ( MN + eb ) - I2 
                - 0.5 * pow(rho,2) * gome_m_sqrd ;

  double alp_3 = MN - ms ;
  double bet_3 = MN * pow(MN - ms, 2) ;
  double gam_3 = pow(MN - ms, 3) ;
  double del_3 = I3 ;
  //..........................................................

  // c = (  (alp_3*bet_1 - alp_1*bet_3) * (alp_2*del_1 - alp_1*del_2) 
  //      - (alp_2*bet_1 - alp_1*bet_2) * (alp_3*del_1 - alp_1*del_3) )
  //     /( (alp_3*bet_1 - alp_1*bet_3) * (alp_2*gam_1 - alp_1*gam_2)
  //       -(alp_2*bet_1 - alp_1*bet_2) * (alp_3*gam_1 - alp_1*gam_3) );
  
  c = (  (alp_3*del_2 - alp_2*del_3) * (alp_2*bet_1 - alp_1*bet_2) 
        - (alp_3*bet_2 - alp_2*bet_3) * (alp_2*del_1 - alp_1*del_2) )
      /( (alp_3*gam_2 - alp_2*gam_3) * (alp_2*bet_1 - bet_2*alp_1)
        -(alp_3*bet_2 - alp_2*bet_3) * (alp_2*gam_1 - alp_1*gam_2) );

  b = ( (alp_2*del_1 - alp_1*del_2) - (alp_2*gam_1 - alp_1*gam_2)*c )
        / (alp_2*bet_1 - alp_1*bet_2) ;

  gsig_m_sqrd = alp_1 / ( del_1 - gam_1*c - bet_1*b) ;

  return *this ;
}
//--------------------------------------------------------------
/// Prints the model parameters
void CompactStar::SigmaOmegaPar::Print() const
{
  std::cout << "\t * ------------------------------------ * \n" ;
  std::cout << "\t |  Sigma - Omega Model Description     | \n" ;
  std::cout << "\t * ------------------------------------ * \n" ;
  char tmp[150] ;
  snprintf(tmp, sizeof(tmp), "\t | %-15s\t %-6.4f\t  fm^2  |\n", "(g_\u03C3/m_\u03C3)^2", gsig_m_sqrd) ;
  std::cout << tmp ;
  snprintf(tmp, sizeof(tmp), "\t | %-15s\t %-6.4f\t  fm^2  |\n", "(g_\u03C9/m_\u03C9)^2", gome_m_sqrd) ;
  std::cout << tmp ;
  snprintf(tmp, sizeof(tmp), "\t | %-15s\t %-6.4f\t  fm^2  |\n", "(g_\u03C1/m_\u03C1)^2", grho_m_sqrd) ;
  std::cout << tmp ;
  snprintf(tmp, sizeof(tmp), "\t | %-15s\t %-6.4f\t        |\n", "b (x100)", b*100) ;
  std::cout << tmp ;
  snprintf(tmp, sizeof(tmp), "\t | %-15s\t %-6.4f\t        |\n", "c (x100)", c*100) ;
  std::cout << tmp ;
  std::cout << "\t * ------------------------------------ * \n" ;

}
//--------------------------------------------------------------

//==============================================================
