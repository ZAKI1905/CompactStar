/*
  Particle Class
*/

#include <Zaki/Math/IntegralTable.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/EOS/Particle.hpp"
#include "CompactStar/EOS/Baryon.hpp"
#include "CompactStar/EOS/Lepton.hpp"


//==============================================================
//             Particle Class
//==============================================================
// Constructor
CompactStar::Particle::Particle(const int& idx, const double& charge, 
                          const double& mass, const double& iso)
  : idx(idx), Q(charge), M(mass), I3(iso)
{
  mtx.lock() ;

  for (auto &&p : P_List)
  {
    if( p->idx == idx)
    {
      Z_LOG_ERROR(("The particle index must be unique, but it is"
                  " the same as '" + p->name + "' index")) ;
      mtx.unlock() ;
      return ;
      // break ;
    }
  }

  P_List.emplace_back(this) ;
  mtx.unlock() ;
  
  // for (auto &&p : L_List)
  // {
  //   if( p->idx == idx)
  //   {
  //     Z_LOG_ERROR(("The lepton index must be unique, but it is"
  //                 " the same as '" + p->name + "' index")) ;
  //     break ;
  //   }
  // }
}
//--------------------------------------------------------------
/// Destructor
CompactStar::Particle::~Particle() 
{
  mtx.lock() ;
  // Removing the pointer from the B_List:
  auto itr = std::find(P_List.begin(), P_List.end(), this) ;
  if(itr != P_List.end())
  {
    P_List.erase(itr) ;
  }
  mtx.unlock() ;
}

//--------------------------------------------------------------
/// Derivative of the chemical potential
double CompactStar::Particle::Dmu(const int& in_idx) 
{
  // Diagonal
  if ( idx != in_idx )
    return 0 ;

  // Assuming gsigma is the first index (not needed really)
  if ( in_idx == 0 )
    return 0 ;
  
  double out = M_PI * M_PI * Baryon::GetBaryonRho() / k() ;
  out       *= 1./ Mu() ;

  return out ;
}
//--------------------------------------------------------------
/// This sets the number density
void CompactStar::Particle::SetRho(const double& rho) 
{
  Rho = rho ;
  number_density_flag = true ;
}

//--------------------------------------------------------------
/// Fermi momentum
double CompactStar::Particle::k() const
{
  return pow(3 * M_PI * M_PI * Rho, 1./3.) ;
}
//--------------------------------------------------------------
/// Fermi momentum squared
double CompactStar::Particle::k2() const
{
  return pow(3 * M_PI * M_PI * Rho, 2./3.) ;
}
//--------------------------------------------------------------
/// Chemical potential
// double CompactStar::Particle::Mu()
// {
//   if (!number_density_flag)
//   {
//     Z_LOG_ERROR("Rho is not set, returning -1.") ;
//     return -1 ;
//   }

//   return sqrt( k2() + M*M ) ;
// }
// Set the chemical potential
void CompactStar::Particle::Set_Mu(const double& mu_in)
{
  mu = mu_in;
}

// Get the chemical potential
double CompactStar::Particle::Mu()
{
    return mu;
}

// Update number density based on chemical potential
void CompactStar::Particle::Set_Rho_From_Mu()
{
    double kF = sqrt(mu * mu - M * M);
    Rho = pow(kF, 3) / (3.0 * M_PI * M_PI);
    number_density_flag = true;
}

// Get number density
double CompactStar::Particle::GetRho() const
{
    return Rho;
}

// // Get particle name
// std::string CompactStar::Particle::Name() const
// {
//     return Get_Name();
// }

//--------------------------------------------------------------
/// Derivative of the chemical potential wrt rho
double CompactStar::Particle::DMu_DRho()
{
  if (!number_density_flag)
  {
    Z_LOG_ERROR("Rho is not set, returning -1.") ;
    return -1 ;
  }

  return M_PI*M_PI / ( Mu() * k() ) ;
}

//--------------------------------------------------------------
/// Contribution to the energy
double CompactStar::Particle::ETerm() const 
{
  return Zaki::Math::I_2({0, k()}, {M}) / M_PI / M_PI ;
}

//--------------------------------------------------------------
/// Contribution to the pressure
double CompactStar::Particle::PTerm() const 
{
  return Zaki::Math::I_3({0, k()}, {M}) / M_PI / M_PI / 3;
}

//--------------------------------------------------------------
/// Prints particle info
void CompactStar::Particle::Print() const 
{
  std::cout << name << ", q= " << Q << ", M= " << M << "\n" ;
}
// //--------------------------------------------------------------
// /// Returns the total baryon density
// double CompactStar::Particle::GetBaryonRho() 
// {
//   double out = 0 ;
//   for (auto &&p : B_List)
//   {
//     if (!p->number_density_flag)
//     {
//       Z_LOG_ERROR("Not all the Rho's are set, returning -1.") ;
//       return -1 ;
//     }

//     out += p->Rho ;
//   }

//   return out ;
// }

//--------------------------------------------------------------
CompactStar::Particle*
CompactStar::Particle::Find(const int& id) 
{
  auto it = std::find_if(P_List.begin(), P_List.end(), 
    [&id](const Particle* obj) {return obj->idx == id;}) ;

  if (it != P_List.end())
    return *it ;
  
  // // Checking the leptons
  // it = std::find_if(L_List.begin(), L_List.end(), 
  //   [&id](const Particle* obj) {return obj->idx == id;}) ;

  // if (it != L_List.end())
  //   return *it ;

  // Not a lepton or baryon
  return nullptr ;
}

//--------------------------------------------------------------
/// Empties the total particle pool (leptons & baryons)
void CompactStar::Particle::EmptyPool() 
{
  Baryon::EmptyBPool() ;
  Lepton::EmptyLPool() ;
}

//--------------------------------------------------------------

//==============================================================
