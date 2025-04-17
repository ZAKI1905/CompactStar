/*
  Sigma-Omega-Rho: npemu Model Class
*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

// Derivative needed for finding 
// the Hesssian of epsilon
#include <gsl/gsl_deriv.h>

#include <Zaki/Math/IntegralTable.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Math/GSLMultiFdfWrapper.hpp>
#include <Zaki/Math/GSLMultiFWrapper.hpp>

#include "CompactStar/ChemicalHeating/SigmaOmegaRho_npemu.hpp"
#include "CompactStar/EOS/Common.hpp"

using namespace Zaki::Physics ;
//==============================================================
//                   Sigma-Omega-Rho Model Class
//==============================================================
// Constructor
CompactStar::SigmaOmegaRho_npemu::SigmaOmegaRho_npemu() 
: Model({0.01, 10*0.153 }),
  proton      (1, +1, PROTON_M_FM, +1./2., 1, 1, 1),
  muon        (2, -1, MUON_M_FM),
  // sigma_minus (3, -1, SIGMA_MINUS_M_FM, -1, 0.568, 0.6, 0.6),
  // lambda      (4,  0, LAMBDA_ZERO_M_FM,  0, 0.568, 0.6, 0.6),
  // sigma_zero  (5,  0, SIGMA_ZERO_M_FM, 0, 0.568, 0.6, 0.6),
  // sigma_plus  (6, +1, SIGMA_PLUS_M_FM, +1, 0.568, 0.6, 0.6),
  // xi_minus    (7, -1, XSI_MINUS_M_FM, -1./2., 0.568, 0.6, 0.6),
  // xi_zero     (8,  0, XSI_ZERO_M_FM, +1./2., 0.568, 0.6, 0.6),
  neutron     (3,  0, NEUTRON_M_FM, -1./2., 1, 1, 1),
  electron    (4, -1, ELECTRON_M_FM)
{
  proton.SetName("Proton") ;
  neutron.SetName("Neutron");
  electron.SetName("Electron") ;
  muon.SetName("Muon") ;
  // sigma_minus.SetName("Sigma-") ;
  // sigma_zero.SetName("Sigma0") ;
  // sigma_plus.SetName("Sigma+") ;
  // xi_minus.SetName("Xi-") ;
  // xi_zero.SetName("Xi0") ;
  // lambda.SetName("Lambda0") ;
  
  // Setting the important rho ranges:
  proton.inclusion      = valid_rho    ;
  neutron.inclusion     = valid_rho    ;
  electron.inclusion    = valid_rho    ;
  muon.inclusion        = {0.13, 1.08} ;
  // sigma_minus.inclusion = {0.27, valid_rho.max} ;
  // lambda.inclusion      = {0.384, valid_rho.max} ;
  // sigma_zero.inclusion  = {0.78, valid_rho.max} ;
  // sigma_plus.inclusion  = {1.05, valid_rho.max} ;
  // xi_minus.inclusion    = {1.35, valid_rho.max} ;
  // // Not sure where xi0 will show up!
  // xi_zero.inclusion     = {1.1*valid_rho.max, 2*valid_rho.max} ;

  SetName("SigmaOmega_npemu") ;
}

//--------------------------------------------------------------
/// Destructor
CompactStar::SigmaOmegaRho_npemu::~SigmaOmegaRho_npemu() 
{
  // std::cout << "-->  " << 1115.683/939 << " =? "<< lambda.M / MN << "\n" ;
}

//--------------------------------------------------------------
void CompactStar::SigmaOmegaRho_npemu::SetPars(const SigmaOmegaPar& in_pars)
{
  params = in_pars ;

  Baryon::SetPars(in_pars) ;
}

//--------------------------------------------------------------
void CompactStar::SigmaOmegaRho_npemu::SetXRhos(const gsl_vector* x)
{

  // Proton density ratio
  // const double xp = gsl_vector_get(x,1);

  // Muon density ratio
  // const double xmu = gsl_vector_get(x,2);
  // const double xmu = 0 ;

  // Sigma- density ratio
  // const double xsig_minus = gsl_vector_get(x,3);
  // const double xsig_minus = 0 ;

  // proton.SetRho(xp * rho) ;
  // sigma_minus.SetRho(xsig_minus * rho) ;
  // neutron.SetRho((1 - xp - xsig_minus) * rho) ;
  // muon.SetRho(xmu * rho) ;
  // electron.SetRho((xp - xmu - xsig_minus ) * rho) ;

  //            Baryons
  for (auto&& b : Baryon::B_List)
  {
    if(b->idx < eq_size)
    {  
      b->SetRho(gsl_vector_get(x, b->idx)*rho) ; 
    }
    else // Baryon number conservation
    {
      double tmp_rho = 0 ;
      for (auto&& b2 : Baryon::B_List)
      {
        if(b2->idx < eq_size)
          tmp_rho += gsl_vector_get(x, b2->idx) ;
      }
      b->SetRho(( 1. - tmp_rho ) *rho) ;  
    }
  }

  //            Leptons
  for (auto&& l : Lepton::L_List)
  {
    if(l->idx < eq_size)
    {  
      l->SetRho(gsl_vector_get(x, l->idx)*rho) ; 
    }
    else // Electric charge conservation
    {
      double tmp_rho = 0 ;
      for (auto&& l2 : Lepton::L_List)
      {
        if(l2->idx < eq_size)
          tmp_rho += l2->Q * gsl_vector_get(x, l2->idx) * rho ;
      }
      // This works assuming that Baryon densities are all set first
      for (auto&& b2 : Baryon::B_List)
      {
        tmp_rho += b2->Q * b2->Rho ;
      }

      l->SetRho( tmp_rho ) ; 
    }
  }
}

//--------------------------------------------------------------
int CompactStar::SigmaOmegaRho_npemu::Eq(const gsl_vector* x, gsl_vector* in_f)
{
  // g_sigma * sigma
  // const double gsig = gsl_vector_get(x,0);

  // // Proton density ratio
  // const double xp = gsl_vector_get(x,1);

  // // Muon density ratio
  // const double xmu = gsl_vector_get(x,2);
  // const double xmu = 0 ;
  // // Sigma- density ratio
  // const double xsig_minus = gsl_vector_get(x,3);
  // const double xsig_minus = 0 ;

  // ...................... SAFTEY CHECK ........................
  // Ratios should be between [ 0 , 1 ]
  // and 'g_sig*sig' should be positive
  for (size_t i = 1; i < eq_size; i++)
  {
    
    if ( gsl_vector_get(x,i) < 0 || gsl_vector_get(x,i) > 1 || gsl_vector_get(x,0) < 0)
    {
      for (size_t j = 0; j < eq_size; j++)
      {
        gsl_vector_set (in_f, j, 10);
      }

      return GSL_SUCCESS ;
    }  
  }
  // ............................................................
  
  // if( xp + xsig_minus > 1 || xp < xsig_minus + xmu )
  // {
  //   for (size_t j = 0; j < eq_size; j++)
  //   {
  //     gsl_vector_set (in_f, j, 10);
  //   }

  //   return GSL_SUCCESS ;
  // }

  // Charge conservation
  // proton.SetRho(xp * rho) ;
  // sigma_minus.SetRho(xsig_minus * rho) ;
  // neutron.SetRho((1 - xp - xsig_minus) * rho) ;
  // muon.SetRho(xmu * rho) ;
  // electron.SetRho((xp - xmu - xsig_minus ) * rho) ;

  SetXRhos(x) ;

  // ...................... SAFTEY CHECK ........................
  // Baryons
  for(auto&& b : Baryon::B_List)
  {
    if( b->Rho < 0 )
    {
      for (size_t j = 0; j < eq_size; j++)
      {
        gsl_vector_set (in_f, j, 11);
      }

      return GSL_SUCCESS ;
    }
  }
  // Leptons
  for(auto&& l : Lepton::L_List)
  {
    if( l->Rho < 0 )
    {
      for (size_t j = 0; j < eq_size; j++)
      {
        gsl_vector_set (in_f, j, 12);
      }
      // std::cout << "Negative rho for '" << l->name << "' .\n" ;
      return GSL_SUCCESS ;
    }
  }
  // ............................................................

  gsig = gsl_vector_get(x, 0) ;

  // double y[eq_size] ;
  
  // // (g_sigma * sig) equation
  // y[0] = Baryon::gsigEq(gsig) ;

  // // Chemical equilibrium equations
  // // Beta decay
  // y[1] = proton.Mu(gsig) - neutron.Mu(gsig) + electron.Mu() ;
 
  // // Muon decay
  // y[2] = muon.Mu() - electron.Mu() ;

  // // N + N -> N + sigma + K
  // y[3] = sigma_minus.Mu(gsig) - neutron.Mu(gsig) - electron.Mu() ;

  SetEquilibriumEqs(rho) ;

  for (size_t i = 0; i < eq_size; i++)
  {
    gsl_vector_set (in_f, i, equilibrium_eqs[i]);
  }

  return GSL_SUCCESS ;
}

//--------------------------------------------------------------
/// Equilibrium Equation Solver
int CompactStar::SigmaOmegaRho_npemu::SolveEq()
{

  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids ;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T, eq_size) ;

  // char tmp_char[100] ;
  // sprintf(tmp_char, "Using '%s' solver...\n", gsl_multiroot_fsolver_name(s)) ;
  // Z_LOG_INFO(tmp_char);


  Zaki::Math::GSLMultiFWrapper<SigmaOmegaRho_npemu, 
    int (SigmaOmegaRho_npemu::*)(const gsl_vector*, gsl_vector*)> 
  func(this, &SigmaOmegaRho_npemu::Eq, eq_size) ;

  gsl_multiroot_function F = static_cast<gsl_multiroot_function> (func) ; 

  // gsl_set_error_handler_off() ;

  gsl_vector *r_guess = gsl_vector_alloc (eq_size);

  FindGuess(rho) ;
  for (size_t i = 0; i < eq_size; i++)
  {
    gsl_vector_set (r_guess, i, root_guess[i]) ;
  }

  gsl_multiroot_fsolver_set (s, &F, r_guess);

  int status;
  int iter = 0, max_iter = 1000;
  // gsl_vector* r ;
  do
  {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);

//     PrintState(iter, s);
    
    if (status)   /* check if solver is stuck */
    {
      Z_LOG_ERROR("Solver is stuck for rho = " + std::to_string(rho) +"!") ;
      break;
    }
    status = gsl_multiroot_test_delta (s->dx, s->x, 1e-10, 1e-10);
    // status = gsl_multiroot_test_residual (s->f, 1e-9) ;
    if(iter == max_iter - 1)
    {
      Z_LOG_ERROR("Solver ran out of iterations, without convergence!") ;
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  // ......................................................
  // Saving the results
  gsig = gsl_vector_get(s->x, 0) ;

  // double xp = gsl_vector_get(s->x, 1) ;
  // double xmu = gsl_vector_get(s->x, 2) ;
  // double xmu = 0 ;
  // double xsig_minus = gsl_vector_get(s->x, 3) ;
  // double xsig_minus = 0 ;

  // proton.SetRho(xp * rho) ;
  // sigma_minus.SetRho(xsig_minus * rho) ;
  // neutron.SetRho((1 - xp - xsig_minus) * rho) ;
  // muon.SetRho(xmu * rho) ;
  // electron.SetRho((xp - xmu - xsig_minus ) * rho) ;

  SetXRhos(s->x) ;
 
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free(r_guess) ;
  
  return status ;
}
//--------------------------------------------------------------
// Prints the status
void CompactStar::SigmaOmegaRho_npemu::PrintState(const size_t& iter, 
                                      gsl_multiroot_fsolver * s)
{
  char tmpc[150] ;
  snprintf(tmpc, sizeof(tmpc), "%2zu) x=", iter) ;

  std::string tmp(tmpc) ;
  for (size_t i = 0; i < eq_size; i++)
  {
    char tmpc2[100] ;
    snprintf(tmpc2, sizeof(tmpc2), "%.1e", gsl_vector_get (s->x, i)) ;
    tmp +=  tmpc2 ;
    if(i != eq_size - 1)
      tmp += ", " ;
  }
  tmp += "\t f=" ;
  for (size_t i = 0; i < eq_size; i++)
  {
    char tmpc2[100] ;
    snprintf(tmpc2, sizeof(tmpc2), "%.1e", gsl_vector_get (s->f, i)) ;
    tmp +=  tmpc2 ;
    if(i != eq_size - 1)
      tmp += ", " ;
  }

  std::cout << tmp << "\n" ;
}

//--------------------------------------------------------------
/// Jacobian
int CompactStar::SigmaOmegaRho_npemu::Jacob(const gsl_vector* x, gsl_matrix* J)
{
  const double gsig = gsl_vector_get(x,0);
  
  // Proton density ratio
  const double xp = gsl_vector_get(x,1);

  // Muon density ratio
  const double xmu = gsl_vector_get(x,2);

  // Sigma- density ratio
  const double xsig_minus = gsl_vector_get(x,3);
  // const double xsig_minus = 0 ;

  // ..................... Saftey .........................
  // Ratios should be between [ 0 , 1 ]
  for (size_t i = 1; i < eq_size; i++)
  {
    
    if ( gsl_vector_get(x,i) < 0 || gsl_vector_get(x,i) > 1)
    {
      for (size_t j = 0; j < eq_size; j++)
      {
        for (size_t k = 0; k < eq_size; k++)
        { 
          gsl_matrix_set (J, k, j, 5)  ;
        }
      }
      return GSL_SUCCESS ;
    }  
  }
  
  if( xp + xsig_minus > 1 || xp < xsig_minus + xmu|| gsl_vector_get(x,0) < 0)
  {
    for (size_t j = 0; j < eq_size; j++)
    {
      for (size_t k = 0; k < eq_size; k++)
      { 
        gsl_matrix_set (J, k, j, 5)  ;
      }
    }
    return GSL_SUCCESS ;
  }
  // ......................................................

  for (auto&& b : Baryon::B_List)
  {
    if(b->idx < eq_size)
    {  
      b->SetRho(gsl_vector_get(x, b->idx)*rho) ; 
    }
    else
    {
      double tmp_rho = 0 ;
      for (auto&& b2 : Baryon::B_List)
      {
        if(b2->idx < eq_size)
          tmp_rho += gsl_vector_get(x, b2->idx) ;
      }
      b->SetRho(( 1. - tmp_rho ) *rho) ; 
    }
  }
  
  for (auto&& l : Lepton::L_List)
  {
    if(l->idx < eq_size)
    {  
      l->SetRho(gsl_vector_get(x, l->idx)*rho) ; 
    }
    else
    {
      double tmp_rho = 0 ;
      for (auto&& l2 : Lepton::L_List)
      {
        if(l2->idx < eq_size)
          tmp_rho += l2->Q * gsl_vector_get(x, l2->idx) * rho ;
      }
      for (auto&& b2 : Baryon::B_List)
      {
        tmp_rho += b2->Q * b2->Rho ;
      }

      l->SetRho( tmp_rho ) ; 
    }
  }

  // Charge conservation
  // proton.SetRho(xp * rho) ;
  // sigma_minus.SetRho(xsig_minus * rho) ;
  // neutron.SetRho((1 - xp - xsig_minus) * rho) ;
  // muon.SetRho(xmu * rho) ;
  // electron.SetRho((xp - xmu - xsig_minus ) * rho) ;

  for (size_t i = 0; i < eq_size; i++)
  {

    // g_sig * sigma equation derivative
    gsl_matrix_set (J, 0, i, Baryon::Dgsig(i, gsig)) ;

    // Beta decay equilibrium condition derivative
    gsl_matrix_set (J, 1, i, 
     proton.Dmu(i, gsig) - neutron.Dmu(i, gsig) + electron.Dmu(i) );
    
    // mu -> e + nu + nu_bar
    // gsl_matrix_set (J, 2, i, 
    //  muon.Dmu(i) - electron.Dmu(i) );

    // n + n -> Sig + K + n
    // gsl_matrix_set (J, 3, i, 
    //  sigma_minus.Dmu(i, gsig) - neutron.Dmu(i, gsig) - electron.Dmu(i) );

  }

  // std::cout << "x = " << xp << ", " << xmu << ", " << xsig_minus << ", " << gsig ;
  // std::cout << "J = " << gsl_matrix_get(J, 0, 0) << " " << gsl_matrix_get(J, 0, 1) << " " << gsl_matrix_get(J, 0, 2) << " " << gsl_matrix_get(J, 0, 3) << "\n" ;
  // std::cout << "J = " << gsl_matrix_get(J, 1, 0) << " " << gsl_matrix_get(J, 1, 1) << " " << gsl_matrix_get(J, 1, 2) << " " << gsl_matrix_get(J, 1, 3) << "\n" ;
  // std::cout << "J = " << gsl_matrix_get(J, 2, 0) << " " << gsl_matrix_get(J, 2, 1) << " " << gsl_matrix_get(J, 2, 2) << " " << gsl_matrix_get(J, 2, 3) << "\n" ;
  // std::cout << "J = " << gsl_matrix_get(J, 3, 0) << " " << gsl_matrix_get(J, 3, 1) << " " << gsl_matrix_get(J, 3, 2) << " " << gsl_matrix_get(J, 3, 3) << "\n" ;
  
  return GSL_SUCCESS ;
}

//--------------------------------------------------------------
int CompactStar::SigmaOmegaRho_npemu::EqFdf(const gsl_vector* x, gsl_vector* in_f, gsl_matrix* J) 
{
  const double gsig = gsl_vector_get(x,0);
  
  // Proton density ratio
  const double xp = gsl_vector_get(x,1);
  
  // Muon density ratio
  const double xmu = gsl_vector_get(x,2);

  // Sigma- density ratio
  const double xsig_minus = gsl_vector_get(x,3);
  // const double xsig_minus = 0 ;

  // ..................... Saftey .........................
  // Ratios should be between [ 0 , 1 ]
  for (size_t i = 1; i < eq_size; i++)
  {
    
    if ( gsl_vector_get(x,i) < 0 || gsl_vector_get(x,i) > 1)
    {
      for (size_t j = 0; j < eq_size; j++)
      {
        for (size_t k = 0; k < eq_size; k++)
        { 
          gsl_matrix_set (J, k, j, 5)  ;
        }
      }

      for (size_t j = 0; j < eq_size; j++)
      { 
        gsl_vector_set (in_f, j, 10);
      }

      return GSL_SUCCESS ;
    }  
  }
   
  if( xp + xsig_minus > 1 || xp < xsig_minus + xmu || gsl_vector_get(x,0) < 0)
  {
    for (size_t j = 0; j < eq_size; j++)
    {
      for (size_t k = 0; k < eq_size; k++)
      { 
        gsl_matrix_set (J, k, j, 5)  ;
      }
    }

    for (size_t j = 0; j < eq_size; j++)
    { 
      gsl_vector_set (in_f, j, 10);
    }

    return GSL_SUCCESS ;
  }
  // ......................................................

  // Charge conservation
  proton.SetRho(xp * rho) ;
  // sigma_minus.SetRho(xsig_minus * rho) ;
  neutron.SetRho((1 - xp - xsig_minus) * rho) ;
  // muon.SetRho(xmu * rho) ;
  electron.SetRho((xp - xmu - xsig_minus ) * rho) ;

  // ........................................................
  //                Evaluating the derivatives
  // ........................................................
  for (size_t i = 0; i < eq_size; i++)
  {

    // g_sig * sigma equation derivative
    gsl_matrix_set (J, 0, i, Baryon::Dgsig(i, gsig)) ;

    // Beta decay equilibrium condition derivative
    gsl_matrix_set (J, 1, i, 
     proton.Dmu(i, gsig) - neutron.Dmu(i, gsig) + electron.Dmu(i) );
    
    // mu -> e + nu + nu_bar
    // gsl_matrix_set (J, 2, i, 
    //  muon.Dmu(i) - electron.Dmu(i) );

    // n + n -> Sig + K + n
    // gsl_matrix_set (J, 3, i, 
    //  sigma_minus.Dmu(i, gsig) - neutron.Dmu(i, gsig) - electron.Dmu(i) );

  }
  // ........................................................

  // ........................................................
  //                Evaluating the function
  // ........................................................
  double y[eq_size] ;
  
  // (g_sigma * sig) equation
  y[0] = Baryon::gsigEq(gsig) ;

  // Checmical equilibrium equations
  // Beta decay
  y[1] = proton.Mu(gsig) - neutron.Mu(gsig) + electron.Mu() ;

  // Muon decay
  // y[2] = muon.Mu() - electron.Mu() ;

  // N + N -> N + sigma + K
  // y[3] = sigma_minus.Mu(gsig) - neutron.Mu(gsig) - electron.Mu() ;



  for (size_t i = 0; i < eq_size; i++)
  {
    gsl_vector_set (in_f, i, y[i]);
  }
  // ........................................................

  // std::cout << "x = " << xp << ", " << xmu << ", " << xsig_minus << ", " << gsig ;
  // std::cout << " y = " << xp << ", " << y[0] << ", " << y[1] << ", " << y[2] << ", " << y[3] << "\n" ;
  // std::cout << "J = " << gsl_matrix_get(J, 0, 0) << " " << gsl_matrix_get(J, 0, 1) << " " << gsl_matrix_get(J, 0, 2) << " " << gsl_matrix_get(J, 0, 3) << "\n" ;
  // std::cout << "J = " << gsl_matrix_get(J, 1, 0) << " " << gsl_matrix_get(J, 1, 1) << " " << gsl_matrix_get(J, 1, 2) << " " << gsl_matrix_get(J, 1, 3) << "\n" ;
  // std::cout << "J = " << gsl_matrix_get(J, 2, 0) << " " << gsl_matrix_get(J, 2, 1) << " " << gsl_matrix_get(J, 2, 2) << " " << gsl_matrix_get(J, 2, 3) << "\n" ;
  // std::cout << "J = " << gsl_matrix_get(J, 3, 0) << " " << gsl_matrix_get(J, 3, 1) << " " << gsl_matrix_get(J, 3, 2) << " " << gsl_matrix_get(J, 3, 3) << "\n" ;
  

  return GSL_SUCCESS ;
}

//--------------------------------------------------------------
/// Equilibrium Equation Solver (using Jacobian)
int CompactStar::SigmaOmegaRho_npemu::SolveEqJacob()
{

  const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridj ;
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc (T, eq_size) ;

  Zaki::Math::GSLMultiFdfWrapper<SigmaOmegaRho_npemu, 
    int (SigmaOmegaRho_npemu::*)(const gsl_vector*, gsl_vector*), 
    int (SigmaOmegaRho_npemu::*)(const gsl_vector*, gsl_matrix*),
    int (SigmaOmegaRho_npemu::*)(const gsl_vector*, gsl_vector*, gsl_matrix*)
    > 
    fdf(this, &SigmaOmegaRho_npemu::Eq, &SigmaOmegaRho_npemu::Jacob, 
      &SigmaOmegaRho_npemu::EqFdf, eq_size) ;


  gsl_multiroot_function_fdf Fdf = static_cast<gsl_multiroot_function_fdf> (fdf) ; 

  // gsl_set_error_handler_off() ;
  gsl_vector *r_guess = gsl_vector_alloc (eq_size);

  FindGuess(rho) ;
  for (size_t i = 0; i < eq_size; i++)
  {
    gsl_vector_set (r_guess, i, root_guess[i]) ;
  }

  gsl_multiroot_fdfsolver_set (s, &Fdf, r_guess);

  int status;
  int iter = 0, max_iter = 1000;
  // gsl_vector* r ;
  do
  {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate (s);

    // ................................................
    // Error handling
    if (status == GSL_EBADFUNC)
    {
      std::cout <<  "Error: The iteration encountered a singular" 
                    " point where the function or its derivative"
                    " evaluated to Inf or NaN.\n" ;
    } else if (status == GSL_EZERODIV)
    {
      std::cout <<  "Error: The derivative of the function vanished"
                    " at the iteration point, preventing the algorithm"
                    " from continuing without a division by zero.\n" ;
    }
    // ................................................

//    PrintState(iter, s);
    
    if (status)   /* check if solver is stuck */
      break;

    // status = gsl_multiroot_test_delta (s->dx, s->x, 1e-10, 1e-10);
    status = gsl_multiroot_test_residual (s->f, 1e-9) ;
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  // ....................................
  // Saving the results
  gsig = gsl_vector_get(s->x, 0) ;

  double xp = gsl_vector_get(s->x, 1) ;
  double xmu = gsl_vector_get(s->x, 2) ;
  double xsig_minus = gsl_vector_get(s->x, 3) ;
  // double xsig_minus = 0 ;


  // proton.SetRho(xp * rho) ;
  // neutron.SetRho((1-xp) * rho) ;
  // electron.SetRho(xp * rho) ;

  proton.SetRho(xp * rho) ;
  // sigma_minus.SetRho(xsig_minus * rho) ;
  neutron.SetRho((1 - xp - xsig_minus) * rho) ;
  // muon.SetRho(xmu * rho) ;
  electron.SetRho((xp - xmu - xsig_minus ) * rho) ;
  // ....................................

  gsl_multiroot_fdfsolver_free(s) ;
  gsl_vector_free(r_guess) ;
  
  return status ;
}
//--------------------------------------------------------------
// Prints the status
void CompactStar::SigmaOmegaRho_npemu::PrintState(const size_t& iter, 
                                      gsl_multiroot_fdfsolver * s)
{
  char tmpc[150] ;
  snprintf(tmpc, sizeof(tmpc), "%2zu) x=", iter) ;

  std::string tmp(tmpc) ;
  for (size_t i = 0; i < eq_size; i++)
  {
    char tmpc2[100] ;
    snprintf(tmpc2, sizeof(tmpc2), "%.1e", gsl_vector_get (s->x, i)) ;
    tmp +=  tmpc2 ;

    if(i != eq_size - 1)
      tmp += ", " ;
  }
  tmp += "\t f=" ;
  for (size_t i = 0; i < eq_size; i++)
  {
    char tmpc2[100] ;
    snprintf(tmpc2, sizeof(tmpc2), "%.1e", gsl_vector_get (s->f, i)) ;
    tmp +=  tmpc2 ;

    if(i != eq_size - 1)
      tmp += ", " ;
  }

  std::cout << tmp << "\n";
}

//--------------------------------------------------------------
// Energy density
double CompactStar::SigmaOmegaRho_npemu::EDens(const double& in_rho) 
{
  SetRho(in_rho) ;
  UpdatePool(in_rho) ;
  SolveEq() ;

  // Sigma cubic term
  double out = params.b*MN*pow(gsig, 3)/3.;

  // Sigma quartic term
  out += params.c*pow(gsig, 4)/4. ;

  // Sigma mass term
  out += pow(gsig,2)/2./params.gsig_m_sqrd ;

  // Omega mass term
  out += 0.5 * pow(Baryon::gpOmega(), 2) / params.gome_m_sqrd ;

  // Rho Isospin term
  out += 0.5 * pow(Baryon::gpRho03(), 2) / params.grho_m_sqrd ;


  // Adding the baryon contributions
  for (auto &&b : Baryon::B_List)
  {
    out += b->ETerm(gsig) ;
  }
  
  // Adding the lepton contributions
  for (auto &&l : Lepton::L_List)
  {
    out += l->ETerm() ;
  }

  out *= INV_FM4_2_G_CM3 ;

  return out ;
}

//--------------------------------------------------------------
// Energy density Hessian
// double CompactStar::SigmaOmegaRho_npemu::EDens(const double& in_rho) 
// {
  // SetRho(in_rho) ;
  // UpdatePool(in_rho) ;
  // SolveEq() ;

  // // Sigma cubic term
  // double out = params.b*MN*pow(gsig, 3)/3.;

  // // Sigma quartic term
  // out += params.c*pow(gsig, 4)/4. ;

  // // Sigma mass term
  // out += pow(gsig,2)/2./params.gsig_m_sqrd ;

  // // Omega mass term
  // out += 0.5 * pow(Baryon::gpOmega(), 2) / params.gome_m_sqrd ;

  // // Rho Isospin term
  // out += 0.5 * pow(Baryon::gpRho03(), 2) / params.grho_m_sqrd ;


  // // Adding the baryon contributions
  // for (auto &&b : Baryon::B_List)
  // {
  //   out += b->ETerm(gsig) ;
  // }
  
//   double out = 0 ;

//   // Adding the lepton contributions
//   for (auto &&l : Lepton::L_List)
//   {
//     out += l->ETerm() ;
//   }

//   out *= INV_FM4_2_G_CM3 ;

//   return out ;
// }

//--------------------------------------------------------------
// Pressure
double CompactStar::SigmaOmegaRho_npemu::Press(const double& in_rho) 
{
  SetRho(in_rho) ;
  UpdatePool(in_rho) ;
  SolveEq() ;

  // Sigma cubic term
  double out = -params.b * MN * pow(gsig, 3)/3.;

  // Sigma quartic term
  out += -params.c * pow(gsig, 4)/4. ;

  // Sigma mass term
  out += -pow(gsig,2) / 2. / params.gsig_m_sqrd ;

  // Omega mass term
  out += 0.5 * pow(Baryon::gpOmega(), 2) / params.gome_m_sqrd ;

  // Rho Isospin term
  out += 0.5 * pow(Baryon::gpRho03(), 2) / params.grho_m_sqrd ;


  // Adding the baryon contributions
  for (auto &&b : Baryon::B_List)
  {
    out += b->PTerm(gsig) ;
  }

  // Adding the lepton contributions
  for (auto &&l : Lepton::L_List)
  {
    out += l->PTerm() ;
  }

  out *= INV_FM4_2_Dyn_CM2 ;

  return out ;
}

//--------------------------------------------------------------
// Updates the pool of leptons & baryons
void CompactStar::SigmaOmegaRho_npemu::UpdatePool(const double& in_rho) 
{
  Particle::EmptyPool() ;

  // neutron.PoolAdd() ;
  // proton.PoolAdd() ;
  // electron.PoolAdd() ;
  // eq_size = 2 ;

  // muon.SetRho(1e-10) ;
  // sigma_minus.SetRho(1e-10) ;
  // lambda.SetRho(1e-10) ;


  // if (muon.inclusion.min < in_rho && in_rho < muon.inclusion.max )
  // {
  //   muon.PoolAdd() ;
  //   eq_size++ ;
  // }

  // if (sigma_minus.inclusion.min < in_rho)
  // {
  //   sigma_minus.PoolAdd() ;
  //   eq_size++ ;
  // }

  // if (lambda.inclusion.min < in_rho)
  // {
  //   lambda.PoolAdd() ;
  //   eq_size++ ;
  // }

  // We are already using the charge & baryon conservation equations
  // so the number of equations is equal to 
  // the number of particles minus 2 plus 1 (+1 is for gsig)
  eq_size = -2 + 1 ;

  for (auto &&p : Particle::P_List)
  {
    // This is just for printing the density of particles not in the pool
    //  [ Could be zero, but in log plots it wouldn't work! ]
    p->SetRho(1e-10) ;

    if (p->inclusion.min <= in_rho && in_rho <= p->inclusion.max )
    {
      p->PoolAdd() ;
      eq_size++ ;
    }
  }
  
  // std::cout << "-> eq_size = " << eq_size << "\n" ;
}

//--------------------------------------------------------------
/// This will generate a row of values for EOS
std::vector<double> CompactStar::SigmaOmegaRho_npemu::EOSRow(const double& rho_i) 
{
  return {EDens(rho_i), Press(rho_i), rho_i, 
          neutron.Rho/rho_i, proton.Rho/rho_i, 
          // sigma_minus.Rho/rho_i,
          // lambda.Rho/rho_i,
          // sigma_zero.Rho/rho_i,
          // sigma_plus.Rho/rho_i,
          // xi_minus.Rho/rho_i,
          // xi_zero.Rho/rho_i,
          electron.Rho/rho_i,
          muon.Rho/rho_i,
    
          neutron.Mu(gsig) / MEV_2_INV_FM,
          proton.Mu(gsig) / MEV_2_INV_FM,
          electron.Mu() / MEV_2_INV_FM,
          muon.Mu() / MEV_2_INV_FM,

          neutron.Dmu(neutron.idx, gsig) / rho_i,
          neutron.Dmu(proton.idx, gsig) / rho_i,
          proton.Dmu(neutron.idx, gsig) / rho_i,
          proton.Dmu(proton.idx, gsig) / rho_i,

          electron.DMu_DRho(),
          muon.DMu_DRho()
          } ;
}

//--------------------------------------------------------------
/// This will generate the header for EOS
std::string CompactStar::SigmaOmegaRho_npemu::EOSHeader() const
{
  char tmp_header[500] ;
  snprintf(tmp_header, sizeof(tmp_header), 
                      "%-14s\t %-14s\t %-14s"
                      "\t %-14s\t %-14s\t %-14s\t %-14s"
                      "\t %-14s\t %-14s\t %-14s\t %-14s"
                      "\t %-14s\t %-14s"
                      "\t %-14s\t %-14s"
                      "\t %-14s\t %-14s"
                      , 
            "e(g/cm^3)", "p(dyne/cm^2)", "rho(1/fm^3)", 
            "r_n", "r_p", "r_e" , "r_mu",
            "Mu_n(MeV)", "Mu_p(MeV)", "Mu_e(MeV)", "Mu_mu(MeV)",
            "d_Mu_nn(fm^2)", "d_Mu_np(fm^2)", 
            "d_Mu_pn(fm^2)", "d_Mu_pp(fm^2)", 
            "d_Mu_ee(fm^2)", "d_Mu_mumu(fm^2)") ;
  return tmp_header ;
}
//--------------------------------------------------------------
/// Contains the information for a good guess depending on rho
void CompactStar::SigmaOmegaRho_npemu::FindGuess(const double& in_rho) 
{

  // In order to avoid bugs rising from removing the muon
  // If rho = rho_1 is such that muon is removed, then for the next
  // rho = rho_2 it will also be removed since the indices are 
  // not reset, if rho_2 > rho_1 this is fine, but if rho_2 < rho_1
  // we might need to include muon again, and rest the indices
  // I'm sure there is a better way to fix this, but for now we
  // just reset all the indices to the initial indices
  proton.idx      =  1 ;
  muon.idx        =  2 ;
  // sigma_minus.idx =  3 ;
  // lambda.idx      =  4 ;
  // sigma_zero.idx  =  5 ;
  // sigma_plus.idx  =  6 ;
  // xi_minus.idx    =  7 ;
  // xi_zero.idx     =  8 ;
  neutron.idx     =  9 ;
  electron.idx    =  10 ;


  if (in_rho < muon.inclusion.min)
    root_guess = { 0.2*MN, 0.03 } ;
  // Muon enters
  else if (in_rho < 0.2)
    root_guess = {0.25*MN, 0.07, 2e-3} ;

  // else if (in_rho < sigma_minus.inclusion.min)
  //   root_guess = {0.3*MN, 0.1, 1e-2} ;

  // // Sigma- enters
  // else if (in_rho < lambda.inclusion.min)
  //   root_guess = {0.45*MN, 0.2, 3e-2, 5e-2} ;

  // // Lambda0 enters
  // else if (in_rho < 0.6)
  //   root_guess = {0.6*MN, 0.23, 1e-3, 0.15, 5e-2} ;
  
  // else if (in_rho < 1.0)
  //   root_guess = {0.6*MN, 0.2, 1e-3, 0.15, 7e-2} ;

  else if (in_rho <= muon.inclusion.max)
    root_guess = {0.715*MN, 0.224, 6e-4} ;
  
  // Muon leaves
  if ( muon.inclusion.max < in_rho ) 
  {
    // Removing muon requires  reordering of the indices:
    // It's lucky that this works! Buggy!!!
    for (auto &&p : Particle::P_List)
    {
      if( p->idx > muon.idx )
      {
        (p->idx)-- ;
      }
    }
    muon.idx = Particle::P_List.size() ;

    // if ( in_rho < sigma_zero.inclusion.min )
    //   root_guess = {0.68*MN, 0.219, 0.214, 3.76e-1} ;
    
    // // Simga0 enters
    // else if ( in_rho < sigma_plus.inclusion.min )
    //   root_guess = {0.68*MN, 0.219, 0.214, 3.76e-1, 1e-2} ;

    // // Simga+ enters
    // else if ( in_rho < xi_minus.inclusion.min )
    //   root_guess = {0.68*MN, 0.219, 0.214, 3.76e-1, 6e-2, 9e-3} ;

    // // Xi- enters
    // else if ( in_rho < xi_zero.inclusion.min )
    //   root_guess = {0.68*MN, 0.219, 0.214, 3.76e-1, 6e-2, 1e-2, 5e-3} ;
    
    // // Xi0 enters
    // else
    //   root_guess = {0.7*MN, 0.21, 0.19, 3e-1, 9e-2, 5e-2, 1e-2, 5e-5} ;
  }


}
//--------------------------------------------------------------
/// Sets all the equations given a specific density 'in_rho' 
/// It sets g_sig*sig equation & chemical equilibrium equations
void CompactStar::SigmaOmegaRho_npemu::SetEquilibriumEqs(const double& in_rho) 
{
  
  equilibrium_eqs.clear() ; 

  // (g_sigma * sig) equation
  equilibrium_eqs.push_back( Baryon::gsigEq(gsig) ) ;

  // .....................................................
  //            Chemical equilibrium equations
  // .....................................................

  // NOTE:
  // This has to change if:
  // 1) we add particles other than leptons or baryons
  // 2) we have baryons with B other than +1, e.g. anti-baryons!
  for (auto &&p : Particle::P_List)
  {
    // Ignore neutron & electron
    if(p->idx == neutron.idx || p->idx == electron.idx )
      continue ;

    if (p->inclusion.min <= in_rho && in_rho <= p->inclusion.max )
    {
      // If it's a baryon:
      Baryon* b = dynamic_cast<Baryon*> (p) ;
      if( b )
        equilibrium_eqs.push_back( b->Mu(gsig) - neutron.Mu(gsig) + b->Q * electron.Mu() ) ;

      // if it's a lepton
      else if (dynamic_cast<Lepton*> (p))
        equilibrium_eqs.push_back( p->Mu() + p->Q * electron.Mu() ) ;

      // otherwise ...
      else
        Z_LOG_ERROR("Particle must be a baryon or lepton!") ;
    }
  }
  

  // Beta decay
  // equilibrium_eqs.push_back( proton.Mu(gsig) - neutron.Mu(gsig) + electron.Mu() ) ;

  // if (muon.inclusion.min <= in_rho && in_rho <= muon.inclusion.max )
  // {
  //   // Muon decay
  //   equilibrium_eqs.push_back( muon.Mu() - electron.Mu()  ) ;
  // }
  // if (sigma_minus.inclusion.min <= in_rho)
  // {
  //   // N + N -> N + sigma + K
  //   equilibrium_eqs.push_back( sigma_minus.Mu(gsig) - neutron.Mu(gsig) - electron.Mu() ) ;
  // }
  
  // if (lambda.inclusion.min <= in_rho)
  // {
  //   equilibrium_eqs.push_back( lambda.Mu(gsig) - neutron.Mu(gsig) ) ;
  // }

  // if (sigma_zero.inclusion.min <= in_rho)
  // {
  //   equilibrium_eqs.push_back( sigma_zero.Mu(gsig) - neutron.Mu(gsig) ) ;
  // }

  // if (sigma_plus.inclusion.min <= in_rho)
  // {
  //   equilibrium_eqs.push_back( sigma_plus.Mu(gsig) - neutron.Mu(gsig) + electron.Mu() ) ;
  // }

  // if (xi_minus.inclusion.min <= in_rho)
  // {
  //   equilibrium_eqs.push_back( xi_minus.Mu(gsig) - neutron.Mu(gsig) - electron.Mu() ) ;
  // }

  // if (xi_zero.inclusion.min <= in_rho)
  // {
  //   equilibrium_eqs.push_back( xi_zero.Mu(gsig) - neutron.Mu(gsig) ) ;
  // }
  // .....................................................

  
}
//--------------------------------------------------------------

//==============================================================
