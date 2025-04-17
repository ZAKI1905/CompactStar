/*
  Lepton Class
*/

#include <Zaki/Math/IntegralTable.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/EOS/Lepton.hpp"

//==============================================================
//             Lepton Class
//==============================================================
// Constructor
CompactStar::Lepton::Lepton(const int& idx, const double& charge,
                            const double& mass)
  : Particle(idx, charge, mass, 0)
{
  // mtx.lock() ;
  // L_List.emplace_back(this) ;
  // mtx.unlock() ;
}
//--------------------------------------------------------------
/// Destructor
CompactStar::Lepton::~Lepton() 
{
  // mtx.lock() ;
  // // Removing the pointer from the L_List:
  // auto itr = std::find(L_List.begin(), L_List.end(), this) ;
  // if(itr != L_List.end())
  // {
  //   L_List.erase(itr) ;
  // }
  // mtx.unlock() ;

  PoolRemove() ;
}
//--------------------------------------------------------------
/// Adds the lepton to the pool of leptons under consideration
/// Useful for dynamically changing the pool without needing
/// heap allocation of particles. 
/// Returns -1 if not successful and 0 if successful
int CompactStar::Lepton::PoolAdd() 
{
  mtx.lock() ;

  // Checking if the lepton already exists:
  auto itr = std::find(L_List.begin(), L_List.end(), this) ;
  if(itr == L_List.end())
  {
    L_List.emplace_back(this) ;
    // Z_LOG_INFO("'"+name+"' was added to the L-pool") ;
    mtx.unlock() ;
    return 0 ;
  }
  else
  {
    // Z_LOG_INFO("'"+name+"' wasn't added to the L-pool"
              //  ", because it's already in it!") ;
    mtx.unlock() ;
    return -1 ;
  }
}

//--------------------------------------------------------------
/// Removes the lepton from the pool of leptons
/// Returns -1 if not successful and 0 if successful
int CompactStar::Lepton::PoolRemove() 
{
  mtx.lock() ;

  // Checking if the lepton already exists:
  auto itr = std::find(L_List.begin(), L_List.end(), this) ;
  if(itr != L_List.end())
  {
    L_List.erase(itr) ;
    // Z_LOG_INFO("'"+name+"' was removed from the L-pool") ;
    mtx.unlock() ;
    return 0 ;
  }
  else
  {
    // Z_LOG_INFO("'"+name+"' wasn't removed from the L-pool"
    //            ", because it's not in it!") ;
    mtx.unlock() ;
    return -1 ;
  }
}

//--------------------------------------------------------------
void CompactStar::Lepton::PrintAll()
{
  for (auto &&l : L_List)
  {
    l->Print() ;
  }
}
//--------------------------------------------------------------
/// Empties the lepton's pool
void CompactStar::Lepton::EmptyLPool() 
{
  L_List.clear() ;
}

//--------------------------------------------------------------

//==============================================================
