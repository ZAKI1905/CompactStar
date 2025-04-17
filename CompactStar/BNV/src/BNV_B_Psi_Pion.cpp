/*
  BNV_B_Psi_Pion class
*/

#include <gsl/gsl_integration.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
// #include <Zaki/Math/Newton.hpp>

#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/BNV/BNV_B_Psi_Pion.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/Core/Pulsar.hpp"

using namespace CompactStar ;

//==============================================================
//                        BNV_B_Psi_Pion class
//==============================================================
// Constructor
BNV_B_Psi_Pion::BNV_B_Psi_Pion(const int& i, const int& j) 
  : 
  BNV_Chi({"B_Psi_Pion_" + std::to_string(i) + std::to_string(j), "B \\to \\psi + \\pi"}),
  gen_i(i), gen_j(j)
  // neutron("neutron", "n", "n", "10", Zaki::Physics::NEUTRON_M_FM,
  //          -3.82608545, 879.6),
  // lambda("lambda", "\\Lambda", "Lam", "100", 
  //         Zaki::Physics::LAMBDA_ZERO_M_FM, -1.22, 2.632e-10), 
  // proton("proton", "p^+", "p", "11", Zaki::Physics::PROTON_M_FM,
  //          5.5856946893, 1e+41)
{ }

//--------------------------------------------------------------
// Destructor
BNV_B_Psi_Pion::~BNV_B_Psi_Pion() { }
//--------------------------------------------------------------
// Returns the time component of Sigma^0 in the baryon CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::get_Sig_0_CM(const BNV_B_Psi_Pion::Params& pars) 
{
  double sig_0_cm = pars.m_B * pars.sig * (pars.x + pars.sig) ;
        sig_0_cm /= sqrt(1 + 2*pars.x*pars.sig + pow(pars.sig, 2)) ;

  return sig_0_cm ;
}

//--------------------------------------------------------------
// Returns the length of 3-vector vec{Sigma} in the baryon CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::get_Sig_V_CM(const BNV_B_Psi_Pion::Params& pars) 
{
  double Sig_V = pars.m_B * pars.sig ;
        Sig_V *= sqrt(-1 + pow(pars.x, 2)) ;
        Sig_V /= sqrt(1 + 2* pars.x * pars.sig + pow(pars.sig,2)) ;

  return Sig_V ;
}

//--------------------------------------------------------------
// Returns the length of 3-vector vec{q} in the baryon CM frame
//    q is the momentum of psi
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::get_q_V_CM(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  // -----------------------------

  double q_v = m_B / (2.*sqrt(1 + 2*x*sig + pow(sig,2))) ;

        q_v *= sqrt( pow(mu_pi,4) + pow(1 - pow(mu_psi,2) + 2*x*sig + pow(sig,2),2) 
                    - 2*pow(mu_pi,2) * (1 + pow(mu_psi,2) + 2*x*sig + pow(sig,2))
                    ) ;
  return q_v ;
}

//--------------------------------------------------------------
// Returns the energy of psi in the baryon CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::get_E_Psi_CM(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  // -----------------------------

  double e_psi  = m_B / 2.0 ;
         e_psi *= 1 - pow(mu_pi,2) + pow(mu_psi,2) + 2*x*sig + pow(sig,2) ;  
         e_psi /= sqrt( 1 + 2*x*sig + pow(sig,2) ) ;

  return e_psi ;
}

//--------------------------------------------------------------
// Returns the energy* of baryon (E_B^*) in the nm frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::get_E_B_s(const BNV_B_Psi_Pion::Params& pars)
{
  return  pars.x * pars.m_B ;
}

//--------------------------------------------------------------
// Returns the energy of baryon (E_B) in the CM frame
// Output checked on Nov 7, 2023 [ failed ]
// Mistake fixed.
double BNV_B_Psi_Pion::get_E_B_CM(const BNV_B_Psi_Pion::Params& pars)
{
  // Wrong Answer:
  // return sqrt( pow(pars.m_B, 2) + 2*get_E_B_s(pars)*pars.sig + pow(pars.sig, 2) ) ;

  // return pars.m_B * sqrt( 1 
  //         + 2*get_E_B_s(pars)*pars.sig/pars.m_B 
  //         + pow(pars.sig, 2) ) ;
   return pars.m_B * sqrt( 1 
          + 2 * pars.x * pars.sig 
          + pow(pars.sig, 2) ) ;         
}

//--------------------------------------------------------------
// Returns the energy* of baryon (E_B^*) in the CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::get_E_B_s_CM(const BNV_B_Psi_Pion::Params& pars)
{
  return pars.m_B * (1 + pars.x*pars.sig) / sqrt(1 + 2*pars.x*pars.sig + pow(pars.sig,2))  ;
}

//--------------------------------------------------------------
//   The angular indep part of the 
//     4-vector dot-product of (q.p*) 
//    evaluated in the CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::q_pstar_zero(const BNV_B_Psi_Pion::Params& pars) 
{
  return get_E_Psi_CM(pars) * get_E_B_s_CM(pars) ;
}

//--------------------------------------------------------------
//   The angular dep part (without cos(theta) ) of the 
//     4-vector dot-product of (q.p*) 
//    evaluated in the CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::q_pstar_first(const BNV_B_Psi_Pion::Params& pars) 
{
  return get_q_V_CM(pars) * get_Sig_V_CM(pars) ;
}

//--------------------------------------------------------------
//   The angular indep part of the 
//     4-vector dot-product of (q*.p*) 
//    evaluated in the CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::qstar_pstar_zero(const BNV_B_Psi_Pion::Params& pars) 
{
  return (get_E_Psi_CM(pars) - get_Sig_0_CM(pars)) * get_E_B_s_CM(pars)  
          - get_Sig_V_CM(pars) *  get_Sig_V_CM(pars) ;
}

//--------------------------------------------------------------
//   The angular indep part of the 
//     4-vector dot-product of (q*.p*) 
//    evaluated in the CM frame
double BNV_B_Psi_Pion::qstar_pstar_first(const BNV_B_Psi_Pion::Params& pars) 
{
  return get_q_V_CM(pars) * get_Sig_V_CM(pars) ;
}

//--------------------------------------------------------------
//   The angular indep part of the 
//     4-vector dot-product of (q.q*) 
//    evaluated in the CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::q_qstar_zero(const BNV_B_Psi_Pion::Params& pars) 
{
  return (get_E_Psi_CM(pars) - get_Sig_0_CM(pars)) * get_E_Psi_CM(pars)  
          - get_q_V_CM(pars) * get_q_V_CM(pars) ;
}

//--------------------------------------------------------------
//   The angular indep part of the 
//     4-vector dot-product of (q.q*) 
//    evaluated in the CM frame
double BNV_B_Psi_Pion::q_qstar_first(const BNV_B_Psi_Pion::Params& pars) 
{
  return get_q_V_CM(pars) * get_Sig_V_CM(pars) ;
}

//--------------------------------------------------------------
//   The angular indep part of the 
//     4-vector dot-product of (q*.q*) 
//    evaluated in the CM frame
// Output verified on Nov 7, 2023
double BNV_B_Psi_Pion::qstar_qstar_zero(const BNV_B_Psi_Pion::Params& pars) 
{
  return pow(get_E_Psi_CM(pars) - get_Sig_0_CM(pars), 2) 
          - get_q_V_CM(pars) * get_q_V_CM(pars) 
          - get_Sig_V_CM(pars) *  get_Sig_V_CM(pars) ;
}

//--------------------------------------------------------------
//   The angular indep part of the 
//     4-vector dot-product of (q*.q*) 
//    evaluated in the CM frame
double BNV_B_Psi_Pion::qstar_qstar_first(const BNV_B_Psi_Pion::Params& pars) 
{
  return 2 * get_q_V_CM(pars) * get_Sig_V_CM(pars) ;
}

//--------------------------------------------------------------
// The W0 term with no extra cos(theta) factors
//   integrated over cos(theta)
double BNV_B_Psi_Pion::phase_W0_zero(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  double T = qstar_qstar_zero(pars) ;
  double q_v = get_q_V_CM(pars) ;
  double Sig_V = get_Sig_V_CM(pars) ;
  // -----------------------------

  double val = -2*b_B*c_B + pow(c_B, 2) ;
  // std::cout << "\t (1) val = " << val << "\n" ;
        
        val += (pow(b_B,2)*(5*pow(m_B,4) - 2*pow(m_B,2)*T + pow(T,2) - 4*pow(q_v,2)*pow(Sig_V,2))) 
                /
                (pow(pow(m_B,2) - T, 2) - 4*pow(q_v,2)*pow(Sig_V,2)) ;
  // std::cout << "\t (2) val = " << val << "\n" ;

  if(q_v * Sig_V > 1e-20)
  {
        val += ( b_B*(b_B - c_B)*pow(m_B, 2) 
                  * log( (pow(m_B,2) - T - 2*q_v*Sig_V)/ (pow(m_B,2) - T + 2*q_v*Sig_V) ) )
                  / ( q_v * Sig_V ) ;
  }
  else
  {
      val += b_B*(b_B - c_B)*pow(m_B, 2) * (-4/(pow(m_B,2) - T));
  }
  // std::cout << "\t (3) Sig_V = " <<  Sig_V  << "\n" ;
  // std::cout << "\t (3) T = " <<  T  << "\n" ;
  // std::cout << "\t (3) val = " <<  val  << "\n" ;

        val *= 2*pow(alpha, 2) / pow(f_pi,2) ;

  // std::cout << "\t (4) val = " << val << "\n" ;


  return val ;
}

//--------------------------------------------------------------
// The W1 term with an extra cos(theta) factor
//   integrated over cos(theta)
double BNV_B_Psi_Pion::phase_W0_first(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  double T = qstar_qstar_zero(pars) ;
  double q_v = get_q_V_CM(pars) ;
  double Sig_V = get_Sig_V_CM(pars) ;
  // -----------------------------

  if(Sig_V < 1e-20)
    return 0 ;
  double first = -4*q_v*(pow(m_B,2) - T)*
                (2*b_B*pow(m_B,2) - c_B*pow(m_B,2) - 
                  b_B*T + c_B*T)*Sig_V ;
  
  double second = 16*(b_B - c_B)*pow(q_v,3)*pow(Sig_V,3) ;

  double third = (-2*b_B*pow(m_B,2) + c_B*pow(m_B,2) + b_B*T - 
                  c_B*T) * (pow(pow(m_B,2) - T,2) - 
                  4*pow(q_v,2)*pow(Sig_V,2))*
                  log((pow(m_B,2) - T - 2*q_v*Sig_V)/
                  (pow(m_B,2) - T + 2*q_v*Sig_V)) ;

  double val = first + second + third ; 
      
      val *= b_B * pow(m_B,2) * pow(alpha, 2) ;

      val /= pow(f_pi,2)*pow(q_v,2)*pow(Sig_V,2)*
              (-pow(pow(m_B,2) - T,2) + 
              4*pow(q_v,2)*pow(Sig_V,2)) ;

  return val ;
}

//--------------------------------------------------------------
// The W1 term with no extra cos(theta) factors
//   integrated over cos(theta)
double BNV_B_Psi_Pion::phase_W1_zero(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  double T = qstar_qstar_zero(pars) ;
  double q_v = get_q_V_CM(pars) ;
  double Sig_V = get_Sig_V_CM(pars) ;
  // -----------------------------

  double val = 8*pow(b_B,2)*pow(m_B,2)*pow(alpha,2) ;
        val /= pow(f_pi,2)*(pow(pow(m_B,2) - T,2) - 
     4*pow(q_v,2)*pow(Sig_V,2)) ;

  return val ;
}

//--------------------------------------------------------------
// The W1 term with an extra cos(theta) factors
//   integrated over cos(theta)
double BNV_B_Psi_Pion::phase_W1_first(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  double T = qstar_qstar_zero(pars) ;
  double q_v = get_q_V_CM(pars) ;
  double Sig_V = get_Sig_V_CM(pars) ;
  // -----------------------------

  double val = 4*q_v*(-pow(m_B,2) + T)*Sig_V ;
// std::cout << "1 val = " << val << "\n" ;
        val += (pow(pow(m_B,2) - T,2) - 
                4*pow(q_v,2)*pow(Sig_V,2))*
                log((-pow(m_B,2) + T - 2*q_v*Sig_V)/
                (-pow(m_B,2) + T + 2*q_v*Sig_V)) ;
// std::cout << "2 val = " << val << "\n" ;

        val *= pow(b_B,2)*pow(m_B,2)*pow(alpha,2) ;
// std::cout << "3 val = " << val << "\n" ;

  if(q_v * Sig_V > 1e-20)
  {
          val /= pow(f_pi,2)*pow(q_v,2)*pow(Sig_V,2)*
                  (-pow(pow(m_B,2) - T,2) + 
                    4*pow(q_v,2)*pow(Sig_V,2)) ;
  }
// std::cout << "4 val = " << val << "\n" ;

  return val ;
}

//--------------------------------------------------------------
// The W1 term with two extra cos(theta) factors
//   integrated over cos(theta)
double BNV_B_Psi_Pion::phase_W1_second(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  double T = qstar_qstar_zero(pars) ;
  double q_v = get_q_V_CM(pars) ;
  double Sig_V = get_Sig_V_CM(pars) ;
  // -----------------------------

  double val = -4*q_v*pow(pow(m_B,2) - T,2)*Sig_V + 
              8*pow(q_v,3)*pow(Sig_V,3) ;

         val += (-pow(m_B,2) + T)*
                (pow(pow(m_B,2) - T,2) - 
                  4*pow(q_v,2)*pow(Sig_V,2))*
                  log((pow(m_B,2) - T - 2*q_v*Sig_V)/
                (pow(m_B,2) - T + 2*q_v*Sig_V)) ;

        val *= pow(b_B,2)*pow(m_B,2)*pow(alpha,2) ;

  if(Sig_V> 1e-20)
  {
        val /= pow(f_pi,2)*pow(q_v,3)*pow(Sig_V,3)*
              (-pow(pow(m_B,2) - T,2) + 
                4*pow(q_v,2)*pow(Sig_V,2)) ;
  }
  return val ;  
}

//--------------------------------------------------------------
// The W0*W1 term with no extra cos(theta) factors
//   integrated over cos(theta)
double BNV_B_Psi_Pion::phase_W0W1_zero(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  double T = qstar_qstar_zero(pars) ;
  double q_v = get_q_V_CM(pars) ;
  double Sig_V = get_Sig_V_CM(pars) ;
  // -----------------------------

  double val = (8*b_B*pow(m_B,2)*q_v*Sig_V)/
                (pow(pow(m_B,2) - T,2) - 
                4*pow(q_v,2)*pow(Sig_V,2)) ;

        val += (-b_B + c_B)*log((-pow(m_B,2) + T - 2*q_v*Sig_V)/
                (-pow(m_B,2) + T + 2*q_v*Sig_V)) ;

        val *= b_B*m_B*pow(alpha,2);
  if(q_v * Sig_V > 1e-20)
  {
        val /= pow(f_pi,2)*q_v*Sig_V ;
  }
  return val ;
}

//--------------------------------------------------------------
// The W0*W1 term with one extra cos(theta) factors
//   integrated over cos(theta)
double BNV_B_Psi_Pion::phase_W0W1_first(const BNV_B_Psi_Pion::Params& pars) 
{
  // -----------------------------
  //      Definitions
  // -----------------------------
  double m_B = pars.m_B ;
  double x = pars.x ;
  double sig = pars.sig ;
  double mu_psi = pars.mu_psi ;
  double mu_pi = pars.mu_pi ;
  double T = qstar_qstar_zero(pars) ;
  double q_v = get_q_V_CM(pars) ;
  double Sig_V = get_Sig_V_CM(pars) ;
  // -----------------------------

  double val = -4*q_v*(pow(m_B,2) - T)*
                (3*b_B*pow(m_B,2) - c_B*pow(m_B,2) - 
                b_B*T + c_B*T)*Sig_V ;

        val += 16*(b_B - c_B)*pow(q_v,3)*pow(Sig_V,3) ;

        val += (3*b_B*pow(m_B,2) - c_B*pow(m_B,2) - b_B*T + 
                c_B*T)*(pow(pow(m_B,2) - T,2) - 
                  4*pow(q_v,2)*pow(Sig_V,2))*
                log((-pow(m_B,2) + T - 2*q_v*Sig_V)/
                (-pow(m_B,2) + T + 2*q_v*Sig_V)) ;

        val *= b_B*m_B*pow(alpha,2) ;
  if(q_v * Sig_V > 1e-20)
  {
        val /= 2*pow(f_pi,2)*pow(q_v,2)*pow(Sig_V,2)*
                (-pow(pow(m_B,2) - T,2) + 
                  4*pow(q_v,2)*pow(Sig_V,2)) ;
  }   
  return val ;
}

//--------------------------------------------------------------
// The A-term in the amplitude
double BNV_B_Psi_Pion::Amp_A(const BNV_B_Psi_Pion::Params& pars) 
{
  // std::cout << "\t phase_W0_zero(pars) = " << phase_W0_zero(pars) << "\n" ;
  // std::cout << "\t q_pstar_zero(pars) = " << q_pstar_zero(pars) << "\n" ;

  // std::cout << "\t phase_W0_first(pars) = " << phase_W0_first(pars) << "\n" ;
  // std::cout << "\t q_pstar_first(pars) = " << q_pstar_first(pars) << "\n" ;

  double A = 2*phase_W0_zero(pars) * q_pstar_zero(pars) ;
        A += 2*phase_W0_first(pars) * q_pstar_first(pars) ;

  double A_W0_0 = 2*phase_W0_zero(pars) * q_pstar_zero(pars) ;

  double A_W0_1 = 2*phase_W0_first(pars) * q_pstar_first(pars) ;

  double A_phase_W0_first = phase_W0_first(pars)  ;
  double A_q_pstar_first = q_pstar_first(pars)  ;
  // ------------------------
  //  Adding W1 factors 
  // ------------------------
  double W1_zero_factor = 2 * qstar_pstar_zero(pars) * q_qstar_zero(pars) 
                         - qstar_qstar_zero(pars) *  q_pstar_zero(pars) ;
  // std::cout << "\t 1 A = " << A << "\n" ;

  double W1_first_factor = 2 * qstar_pstar_zero(pars) * q_qstar_first(pars) 
                          + 2 * qstar_pstar_first(pars) * q_qstar_zero(pars) 
                         - qstar_qstar_zero(pars) *  q_pstar_first(pars) 
                         - qstar_qstar_first(pars) *  q_pstar_zero(pars) ;
  // std::cout << "\t 2 A = " << A << "\n" ;

  double W1_second_factor = 2 * qstar_pstar_first(pars) * q_qstar_first(pars) 
                         - qstar_qstar_first(pars) *  q_pstar_first(pars)  ;

        A += 2*phase_W1_zero(pars) * W1_zero_factor ;

        // std::cout << "\t 3 A = " << A << "\n" ;
        A += 2*phase_W1_first(pars) * W1_first_factor ;
        // std::cout << "\t 4 A = " << A << "\n" ;
        A += 2*phase_W1_second(pars) * W1_second_factor ;

  double A_W1_0 = 2*phase_W1_zero(pars) * W1_zero_factor ;
  double A_W1_1 = 2*phase_W1_first(pars) * W1_first_factor ;
  double A_W1_2 = 2*phase_W1_second(pars) * W1_second_factor ;
  // ------------------------
  // std::cout << "\t 5 A = " << A << "\n" ;


  // ------------------------
  //  Adding W0*W1 factors 
  // ------------------------
        A += -4*phase_W0W1_zero(pars) * pars.m_B * q_qstar_zero(pars) ;
        A += -4*phase_W0W1_first(pars) * pars.m_B * q_qstar_first(pars) ;

  double A_W0W1_0 = -4*phase_W0W1_zero(pars) * pars.m_B * q_qstar_zero(pars) ;
  double A_W0W1_1 = -4*phase_W0W1_first(pars) * pars.m_B * q_qstar_first(pars) ;
  // ------------------------
  // std::cout << "\t 6 A = " << A << "\n" ;

  if (A < 0)
  {
    std::cout << "\n *****  Warning!! ***** \n" ;
    std::cout << "\n A = " << A << "\n" ;
    std::cout << " pars = {mB=" << pars.m_B << ", mu_pi="
              << pars.mu_pi << ", mu_psi="<< pars.mu_psi 
              << ", sig="<< pars.sig << ", x=" << pars.x << "}\n" ;
    std::cout << " A_W0_0=" << A_W0_0 << ", A_phase_W0_first=" << A_phase_W0_first 
              << ", A_q_pstar_first=" << A_q_pstar_first << ", A_W0_1=" << A_W0_1
              << ", A_W1_0=" << A_W1_0 << ", A_W1_1=" << A_W1_1 
              << ", A_W1_2=" << A_W1_2 << ", A_W0W1_0="<< A_W0W1_0 
              << ", A_W0W1_1="<< A_W0W1_1 << "\n" ;
    std::cout << "\n *****====================***** \n" ;
  }

  return A ;
}

//--------------------------------------------------------------
//  Loop integral in dimension-less form
double BNV_B_Psi_Pion::LoopIntegral(const double& xD,const double& xU, const double& xY) 
{
  double val = (4 - xD) * xD * log(xD) / ((1 - xD) * (xD - xU) * (xD - xY));
        val += (4 - xU) * xU * log(xU) / ((1 - xU) * (xU - xD) * (xU - xY));
        val += (4 - xY) * xY * log(xY) / ((1 - xY) * (xY - xD) * (xY - xU));
 
    return val ;
}

//--------------------------------------------------------------
// The averaged amplitude sqrd (integrated over Cos[theta])
double BNV_B_Psi_Pion::Amplitude(const BNV_B_Psi_Pion::Params& pars) 
{ 
  // i and j start from 1 and go to 3.
  size_t i = gen_i, j = gen_j ;

  double amp = GF ;
        amp /= 16. * pow(M_PI*mW, 2) ;
        amp *= V_CKM[i-1][0] * V_CKM[0][j-1] ;
        amp *= m_D[j-1] * m_U[i-1] ;
        amp *= LoopIntegral(pow(m_D[j-1]/mW, 2),
                            pow(m_U[i-1]/mW, 2),
                            pow(m_X/mW, 2)) ;

        amp *= amp ;

        amp *= Amp_A(pars) ;
  
  return amp ;
}

//--------------------------------------------------------------
// The baryon decay rate in the CM frame (in MeV)
double BNV_B_Psi_Pion::CM_Decay_Rate(const BNV_B_Psi_Pion::Params& pars) 
{
  double cm_rate = Amplitude(pars) ;

        cm_rate *= get_q_V_CM(pars) ;

        cm_rate /= get_E_B_CM(pars) * get_E_B_s_CM(pars) ;

        cm_rate /= 16. * M_PI ;

  return cm_rate ;
}

//--------------------------------------------------------------
// The baryon decay rate in the CM frame (in MeV)
double BNV_B_Psi_Pion::NM_Decay_Rate(const BNV_B_Psi_Pion::Params& pars) 
{
  double nm_rate = CM_Decay_Rate(pars) ;
    
        // boost factor
        nm_rate *= get_E_B_s_CM(pars) / get_E_B_s(pars) ;

  return nm_rate ;
}

//--------------------------------------------------------------
// Baryon density loss rate (the output is positive)
double BNV_B_Psi_Pion::PhaseSpace_Integral(const double& x_F) 
{
  // selected_pars = pars ;
  // selected_pars.SetPars(636.774, 0.219858, 0.785208, 3.76196) ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
  double err, result;

  Zaki::Math::GSLFuncWrapper<BNV_B_Psi_Pion, double (BNV_B_Psi_Pion::*)(const double&)> 
    Fp(this, &BNV_B_Psi_Pion::PhaseSpace_Integrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
  gsl_integration_workspace_free(w);

  return result ;
}

//--------------------------------------------------------------
// Phase space integrand function
double BNV_B_Psi_Pion::PhaseSpace_Integrand(const double& x) 
{

  // selected_pars.SetPars(636.774, 0.219858, 0.785208, 3.76196) ;
  selected_pars.x = x ;
  
  if ( get_E_B_CM(selected_pars) <= 
          (selected_pars.mu_psi + selected_pars.mu_pi)*selected_pars.m_B 
          ) 
  {
    return 0 ;
  }

  double n_dot_p = NM_Decay_Rate(selected_pars) ;
  // if (n_dot_p < 0)
  // {
  //   std::cout << "\n *****  Warning!! ***** \n" ;
  //   std::cout << "\n n_dot_p = " << n_dot_p << "\n" ;
  //   std::cout << " pars = {mB=" << selected_pars.m_B << ", mu_pi="
  //             << selected_pars.mu_pi << ", mu_psi="<< selected_pars.mu_psi 
  //             << ", sig="<< selected_pars.sig << ", x=" << selected_pars.x << "}\n" ;
  //   std::cout << "\n *****====================***** \n" ;
  // }
  // factor of 2 for baryon spin, 
  //  and a factor of 2 for the fact that 
  //    subsequent to this decay another baryon 
  //    will be annihilated too, so Delta(B) = -2
  n_dot_p *= 4 ;

  n_dot_p /= 2*M_PI*M_PI ;

  n_dot_p *= pow(selected_pars.m_B, 3) 
              * sqrt(x*x -1) 
              * x ;

  return n_dot_p ;
}

//--------------------------------------------------------------
/// The decay rate per unit volume
/// Out put is in s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Psi_Pion::Rate_vs_Density(
                          const double& m_psi, 
                          const Baryon& B, 
                          const bool& gen_plots)
{  
  Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI*n_B[B.label]).pow(2.0/3.0) ;

  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  Zaki::Vector::DataColumn m_B = m_B_ds[B.label] ;
  Zaki::Vector::DataColumn Sigma_0 = sig0_ds[B.label] ;
  Zaki::Vector::DataColumn x_F = (m_B*m_B + kF_2).sqrt() / m_B ;

  Zaki::Vector::DataSet rate ;
  rate.Reserve(2, n_B["n_tot"].Size() ) ;

  for (size_t i = 0; i < m_B.Size(); i++)
  {
    selected_pars.SetPars(m_B[i], m_pion/m_B[i],  m_psi/m_B[i], Sigma_0[i]/m_B[i] ) ;
    // std::cout << "\t i=" << i << ") pars = {" << selected_pars.m_B << ", "
    //           << selected_pars.mu_pi << ", "<< selected_pars.mu_psi 
    //           << ", "<< selected_pars.sig << ", xF = " << x_F[i] << "}\n" << std::flush ;
    rate.AppendRow({n_B["n_tot"][i], PhaseSpace_Integral(x_F[i]) }) ;
  }
  
  // Converting "MeV^4" into "MeV/fm^3"
  rate[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "MeV/fm^3" into "s^-1/fm^3"
  rate[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  if(gen_plots)
  {
    hidden_Plot_Rate_vs_Density(B, m_psi, &rate) ;
  }

  return rate ;
}

//--------------------------------------------------------------
// Returns rate 
//  as a function of radius in units of s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Psi_Pion::Rate_vs_R(const double& m_psi, 
                                const Baryon& B, 
                                const bool& gen_plots) 
{
  Zaki::Vector::DataSet rate_vs_n = Rate_vs_Density(m_psi, B) ;
  rate_vs_n.Interpolate(0, 1) ;

  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;
  
  Zaki::Vector::DataSet rate_vs_r({pulsar.GetProfile()->operator[](0),
                                   rate_vs_n.Evaluate(1, n_r)}) ;


  if (gen_plots)
  {
    hidden_Plot_Rate_vs_R(B, m_psi, &rate_vs_r) ;
  }
  
  return rate_vs_r ;
}

// //--------------------------------------------------------------
// // Returns rate 
// //  as a function of radius in units of s^-1/fm^3
// Zaki::Vector::DataSet BNV_B_Chi_Photon::Rate_vs_R_Slow(const double& m_chi, 
//                                 const Baryon& B, 
//                                 const bool& gen_plots) 
// {
//   Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI
//                         * micro_r["n_" + B.label]).pow(2.0/3.0) ;
//   kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

//   Zaki::Vector::DataColumn m_B = micro_r["m_" + B.label] ;
//   Zaki::Vector::DataColumn Sigma_0 = micro_r["sig_" + B.label]  ;
//   Zaki::Vector::DataColumn x_F = (m_B*m_B + kF_2).sqrt() / m_B ;

//   Zaki::Vector::DataSet rate_vs_r ;
//   rate_vs_r.Reserve(2, m_B.Size() ) ;

//   for (size_t i = 0; i < m_B.Size(); i++)
//   {
//     rate_vs_r.AppendRow({micro_r[0][i], 
//       PhaseSpace_Integral(B, m_B[i], m_chi, Sigma_0[i], x_F[i]) }) ;
//   }
  
//   // Converting "MeV^4" into "MeV/fm^3"
//   rate_vs_r[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

//   // Converting "MeV/fm^3" into "s^-1/fm^3"
//   rate_vs_r[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

//   if (gen_plots)
//   {
//     hidden_Plot_Rate_vs_R(B, m_chi, &rate_vs_r) ;
//   }
  
//   return rate_vs_r ;
// }

// //--------------------------------------------------------------

//==============================================================   