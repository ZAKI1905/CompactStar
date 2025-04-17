/*
  Baryon Class
*/

// #include <gsl/gsl_math.h>
#include <Zaki/Math/IntegralTable.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/EOS/Baryon.hpp"
#include "CompactStar/EOS/Lepton.hpp"

//==============================================================
//             Baryon Class
//==============================================================
// Constructor
CompactStar::Baryon::Baryon(const int& idx, const double& charge, 
  const double& mass, const double& isospin)
  : Particle(idx, charge, mass, isospin)
{
  // mtx.lock() ;
  // B_List.emplace_back(this) ;
  // mtx.unlock() ;
}
//==============================================================
// Constructor
CompactStar::Baryon::Baryon(const int& idx, const double& charge, 
  const double& mass, const double& isospin, const double& in_x_ome,
  const double& in_x_rho, const double& in_x_sig)
  : Particle(idx, charge, mass, isospin), 
   x_ome(in_x_ome), x_rho(in_x_rho), x_sig(in_x_sig)
{
  // mtx.lock() ;
  // B_List.emplace_back(this) ;
  // mtx.unlock() ;
}
//--------------------------------------------------------------
/// Destructor
CompactStar::Baryon::~Baryon() 
{
  // mtx.lock() ;
  // // Removing the pointer from the B_List:
  // auto itr = std::find(B_List.begin(), B_List.end(), this) ;
  // if(itr != B_List.end())
  // {
  //   B_List.erase(itr) ;
  // }
  // mtx.unlock() ;

  PoolRemove() ;
} 

//--------------------------------------------------------------
/// Adds the baryon to the pool of baryons under consideration
/// Useful for dynamically changing the pool without needing
/// heap allocation of particles. 
/// Returns -1 if not successful and 0 if successful
int CompactStar::Baryon::PoolAdd() 
{
  mtx.lock() ;

  // Checking if the baryon already exists:
  auto itr = std::find(B_List.begin(), B_List.end(), this) ;
  if(itr == B_List.end())
  {
    B_List.emplace_back(this) ;
    // Z_LOG_INFO("'"+name+"' was added to the B-pool") ;
    mtx.unlock() ;
    return 0 ;
  }
  else
  {
    // Z_LOG_INFO("'"+name+"' wasn't added to the B-pool"
    //            ", because it's already in it!") ;
    mtx.unlock() ;
    return -1 ;
  }
}

//--------------------------------------------------------------
/// Removes the Baryon from the pool of Baryons
/// Returns -1 if not successful and 0 if successful
int CompactStar::Baryon::PoolRemove() 
{
  mtx.lock() ;

  // Checking if the baryon already exists:
  auto itr = std::find(B_List.begin(), B_List.end(), this) ;
  if(itr != B_List.end())
  {
    B_List.erase(itr) ;
    // Z_LOG_INFO("'"+name+"' was removed from the B-pool") ;
    mtx.unlock() ;
    return 0 ;
  }
  else
  {
    // Z_LOG_INFO("'"+name+"' wasn't removed from the B-pool"
    //            ", because it's not in it!") ;
    mtx.unlock() ;
    return -1 ;
  }
}

//--------------------------------------------------------------
/// Setting the overall parameters of the model
void CompactStar::Baryon::SetPars(const SigmaOmegaPar& in_pars)
{
  pars = in_pars ;
}

//--------------------------------------------------------------
// /// This sets the B-number density
// void CompactStar::Baryon::SetRho(const double& rho_B) 
// {
//   RhoB = rho_B ;
//   number_density_flag = true ;
// }
//--------------------------------------------------------------
/// Returns the value of [ g_omega(proton) * omega_0 ]
double CompactStar::Baryon::gpOmega() 
{
  double out = 0 ;
  for (auto &&b : B_List)
  {
    // Baryon* b = dynamic_cast<Baryon*> (p) ;
    if (!b->number_density_flag)
    {
      Z_LOG_ERROR("Not all the Rho's are set, returning -1.") ;
      return -1 ;
    }

    out += b->x_ome * b->Rho ;
  }

  return out * pars.gome_m_sqrd ;
}

//--------------------------------------------------------------
/// Returns the value of [ g_rho(proton) * rho_03 ]
double CompactStar::Baryon::gpRho03() 
{
  double out = 0 ;
  for (auto &&b : B_List)
  {
    // Baryon* b = dynamic_cast<Baryon*> (p) ;
    if (!b->number_density_flag)
    {
      Z_LOG_ERROR("Not all the Rho's are set, returning -1.") ;
      return -1 ;
    }

    out += b->x_rho * b->I3 * b->Rho ;
  }

  return out * pars.grho_m_sqrd ;
}

//--------------------------------------------------------------
/// Returns the total baryon density
double CompactStar::Baryon::GetBaryonRho() 
{
  double out = 0 ;
  for (auto &&p : B_List)
  {
    if (!p->number_density_flag)
    {
      Z_LOG_ERROR("Not all the Rho's are set, returning -1.") ;
      return -1 ;
    }

    out += p->Rho ;
  }

  return out ;
}

//--------------------------------------------------------------
/// Chemical potential
double CompactStar::Baryon::Mu(const double& gsig) const
{
  if (!number_density_flag)
  {
    Z_LOG_ERROR("Not all the Rho's are set, returning -1.") ;
    return -1 ;
  }
  
  double out = 0 ;

  out += x_ome * gpOmega() + x_rho * gpRho03() * I3 ;
  out += sqrt( k2() + pow(M - x_sig * gsig, 2) );

  return out ;
}

//--------------------------------------------------------------
/// Equation describing g_sigma * sigma
double CompactStar::Baryon::gsigEq(const double& gsig)
{

  using namespace Zaki::Math ;

  double out = 0 ;

  for (auto &&b : B_List)
  {
    // Baryon* b = dynamic_cast<Baryon*> (p) ;
    if (!b->number_density_flag)
    {
      Z_LOG_ERROR("Not all the Rho's are set, returning -1.") ;
      return -1 ;
    }

    out  += b->x_sig * ( b->M - b->x_sig * gsig ) 
            * I_1({0, b->k()}, {b->M - b->x_sig * gsig}) ;
  }

  out *= 1./M_PI/M_PI ;

  out +=  - pars.b * MN * pow(gsig, 2) ;
  out +=  - pars.c * pow(gsig, 3) ;

  out += - gsig / pars.gsig_m_sqrd  ;

  return out ;
}

//--------------------------------------------------------------
void CompactStar::Baryon::PrintAll()
{
  for (auto &&b : B_List)
  {
    b->Print() ;
  }
}
//--------------------------------------------------------------
// /// Prints baryon info
// void CompactStar::Baryon::Print() const 
// {
//   std::cout << name << ", q= " << Q << ", M= " << M << "\n" ;
// }
//--------------------------------------------------------------
/// Contribution to the energy
double CompactStar::Baryon::ETerm(const double& gsig) const 
{
  return Zaki::Math::I_2({0, k()}, {M - x_sig * gsig}) / M_PI / M_PI ;
}

//--------------------------------------------------------------
/// Contribution to the pressure
double CompactStar::Baryon::PTerm(const double& gsig) const 
{
  return Zaki::Math::I_3({0, k()}, {M - x_sig * gsig}) / M_PI / M_PI / 3;
}

//--------------------------------------------------------------
/// Derivative of the chemical potential
double CompactStar::Baryon::Dmu(const int& in_idx, 
                                const double& gsig) 
{
  Particle* p_i = Find(in_idx) ;

  // Leptons are decoupled
  if ( dynamic_cast<Lepton*> (p_i) )
    return 0 ;

  // Assuming gsigma is the first index:
  if ( in_idx == 0 )
  {
    double tmp = (M - x_sig * gsig ) * x_sig ;
          tmp *= 1. / sqrt( k2() + pow(M - x_sig * gsig, 2) ) ;
    return tmp ;
  }

  //..............................................
  //          It is a baryon index 
  //..............................................
  Baryon* b_i = dynamic_cast<Baryon*> (p_i) ;

  double d_mu_d_x = M_PI * M_PI * GetBaryonRho() / k() ;
  d_mu_d_x *= 1./ sqrt( k2() + pow(M - x_sig * gsig, 2) ) ;

  double d_mu_d_omega = x_ome * pars.gome_m_sqrd ;
        d_mu_d_omega *= b_i->x_ome * b_i->Rho ;

  double d_mu_d_rho = x_rho * I3 * pars.grho_m_sqrd ;
        d_mu_d_rho *= b_i->x_rho * b_i->I3 * b_i->Rho ;

  double out = d_mu_d_x + d_mu_d_omega + d_mu_d_rho ;

  return out ;
}
//--------------------------------------------------------------
/// Derivative of the g_sig * sigma equation
double CompactStar::Baryon::Dgsig(const int& in_idx, 
                                  const double& gsig) 
{
  using namespace Zaki::Math ;
  Particle* p_i = Find(in_idx) ;

  // Leptons are decoupled
  if ( dynamic_cast<Lepton*> (p_i) )
    return 0 ;

  // Assuming gsigma is the first index:
  if ( in_idx == 0 )
  {
    double tmp_sum = 0 ;

    for (auto &&p : B_List)
    {
      Baryon* b = dynamic_cast<Baryon*> (p) ;
      double mb_xsig = b->M - b->x_sig*gsig ;

      tmp_sum += -pow(b->x_sig, 2)
                  * I_1({0, b->k()}, {mb_xsig}) ;
      tmp_sum += -pow(b->x_sig, 2) * mb_xsig
                  * D_I_1_a(b->k(), mb_xsig) ;
    }
    

    double tmp = -2*pars.b * MN *gsig - 3*pars.c * pow(gsig, 2) ;
          tmp += -1./pars.gsig_m_sqrd ;
          tmp += tmp_sum / pow(M_PI,2) ;

    return tmp ;
  }

  //..............................................
  //          It is a baryon index 
  //..............................................
  Baryon* b_i = dynamic_cast<Baryon*> (p_i) ;

 double out = I_1( {0, b_i->k()}, { b_i->M - b_i->x_sig * gsig } ) ;
        out *= b_i->M  -  2 * b_i->x_sig * gsig ;
        out *= 1. / M_PI / M_PI ;

        out += (b_i->M - b_i->x_sig * gsig) * b_i->x_sig * Baryon::GetBaryonRho()
                / sqrt(b_i->k2() + pow(b_i->M - b_i->x_sig * gsig,2) ) ;

  return out ;
}
//--------------------------------------------------------------
double CompactStar::Baryon::D_I_1_a(const double& k, const double& a) 
{
  double out = a * k / sqrt(a*a + k*k) ;
  out       += -a * atanh( k / sqrt(a*a + k*k)) ;

  return out ;
}

//--------------------------------------------------------------
/// Empties the baryon's pool
void CompactStar::Baryon::EmptyBPool() 
{
  B_List.clear() ;
}

//--------------------------------------------------------------

//==============================================================
