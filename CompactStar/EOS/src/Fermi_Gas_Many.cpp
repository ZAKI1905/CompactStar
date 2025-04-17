/*
  Fermi_Gas_Many Model Class Implementation
  Author: Mohammadreza Zakeri (Zaki)
  Date: Nov 23, 2024
  Description:
    Implements the methods declared in Fermi_Gas_Many.hpp to model
    the EOS of neutron star matter with multiple fermion species.
*/

#include "CompactStar/EOS/Fermi_Gas_Many.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Zaki/Physics/Constants.hpp>

using namespace CompactStar;
using namespace Zaki::Physics;

//--------------------------------------------------------------
// Constructor
Fermi_Gas_Many::Fermi_Gas_Many()
    : Model({1e-12, 10e+0}),
      proton(1, +1, PROTON_M_MEV, +0.5),
      neutron(2, 0, NEUTRON_M_MEV, -0.5),
      electron(3, -1, ELECTRON_M_MEV),
      // muon(4, -1, MUON_M_MEV),
      eq_size(1)
{
    proton.SetName("Proton");
    neutron.SetName("Neutron");
    electron.SetName("Electron");
    // muon.SetName("Muon");

    // Add particles to the pools
    baryon_pool = {&proton, &neutron};
    // lepton_pool = {&electron, &muon};
    lepton_pool = {&electron};

    // Initialize root guess and equilibrium equations
    root_guess.resize(eq_size, 0.0);
    equilibrium_eqs.resize(eq_size, 0.0);

    SetName("Fermi_Gas_Many");
}

//--------------------------------------------------------------
// Destructor
Fermi_Gas_Many::~Fermi_Gas_Many() {}

//--------------------------------------------------------------
// Update particle properties based on chemical potentials
void Fermi_Gas_Many::UpdateParticles()
{
  // Update baryons
  for (auto& baryon : baryon_pool)
  {
    baryon->Set_Rho_From_Mu();
  }

  // Update leptons
  for (auto& lepton : lepton_pool)
  {
    lepton->Set_Rho_From_Mu();
  }
}

//--------------------------------------------------------------
// Set up equilibrium equations
void Fermi_Gas_Many::SetEquilibriumEqs()
{
  double rho_n = root_guess[0]; // rho_n
  std::cout << "rho = " << rho ;
  chem_equilibrium_eq = sqrt(pow(neutron.M * Zaki::Physics::MEV_2_INV_FM, 2) + pow(3*M_PI*M_PI*rho_n, 2/3)) 
                        - sqrt(pow(proton.M* Zaki::Physics::MEV_2_INV_FM, 2) + pow(3*M_PI*M_PI*(rho - rho_n), 2/3)) 
                        - sqrt(pow(electron.M* Zaki::Physics::MEV_2_INV_FM, 2) + pow(3*M_PI*M_PI*(rho - rho_n), 2/3)) ;

  // // Extract chemical potentials from root_guess
  // double mu_e = root_guess[0];
  // double mu_p = root_guess[1];

  // std::cout << " mu_e = "  << mu_e << "\n" ;
  // std::cout << " mu_p = "  << mu_p << "\n" ;

  // // Set chemical potentials
  // electron.Set_Mu(mu_e);
  // // muon.Set_Mu(mu_e); // Muon chemical potential equals electron's
  // proton.Set_Mu(mu_p);
  // neutron.Set_Mu(mu_p + mu_e); // From beta equilibrium

  // // Update particle densities
  // UpdateParticles();

  // // Equations:
  // // 1. Charge neutrality: Q_total = 0
  // equilibrium_eqs[0] = proton.GetRho() * proton.Q +
  //                       neutron.GetRho() * neutron.Q +
  //                       electron.GetRho() * electron.Q ;
  //                     //  + muon.GetRho() * muon.Q;

  // // 2. Baryon number conservation: rho = rho_p + rho_n
  // equilibrium_eqs[1] = proton.GetRho() + neutron.GetRho() - rho;

  std::cout << "chem_equilibrium_eq = " << chem_equilibrium_eq << "\n" ;
  // std::cout << "equilibrium_eqs[1] = " << equilibrium_eqs[1] ;
}

//--------------------------------------------------------------
// Provide an initial guess for the solver
void Fermi_Gas_Many::FindGuess()
{
  // Simple initial guess based on non-interacting Fermi gas
  // root_guess[0] = electron.M ; // Electron chemical potential
  // root_guess[1] = proton.M ;   // Proton chemical potential
  // root_guess[0] = proton.M ; // Electron chemical potential
  // root_guess[1] = proton.M ;   // Proton chemical potential
  root_guess[0] = 0.7* rho ;
}

//--------------------------------------------------------------
// GSL wrapper for the equilibrium function
int Fermi_Gas_Many::EquilibriumFunction(const gsl_vector* x, void* params, gsl_vector* f)
{
  // Fermi_Gas_Many* model = static_cast<Fermi_Gas_Many*>(params);

  // // Extract variables
  // model->root_guess[0] = gsl_vector_get(x, 0); // mu_e
  // model->root_guess[1] = gsl_vector_get(x, 1); // mu_p

  // // Set equilibrium equations
  // model->SetEquilibriumEqs();

  // // Set function values
  // gsl_vector_set(f, 0, model->equilibrium_eqs[0]);
  // gsl_vector_set(f, 1, model->equilibrium_eqs[1]);

  // return GSL_SUCCESS;

  Fermi_Gas_Many* model = static_cast<Fermi_Gas_Many*>(params);

  // Extract variables
  model->root_guess[0] = gsl_vector_get(x, 0); // rho_n
  // model->root_guess[1] = gsl_vector_get(x, 1); // mu_p

  // Set equilibrium equations
  model->SetEquilibriumEqs();

  // double rho_n = gsl_vector_get(x, 0); // rho_n
  // double chem_equilibrium_eq = sqrt(neutron.M ** 2 + pow(3*M_PI*M_PI*rho_n, 2/3)) 
  //                       - sqrt(proton.M ** 2 + pow(3*M_PI*M_PI*(rho - rho_n), 2/3)) 
  //                       - sqrt(electron.M ** 2 + pow(3*M_PI*M_PI*(rho - rho_n), 2/3)) ;

  // Set function values
  gsl_vector_set(f, 0, model->chem_equilibrium_eq);
  // gsl_vector_set(f, 1, model->equilibrium_eqs[1]);

  return GSL_SUCCESS;
}

//--------------------------------------------------------------
// // Solve the equilibrium equations
// int Fermi_Gas_Many::SolveEquations()
// {
//     const size_t n = eq_size;

//     gsl_multiroot_function f = {&EquilibriumFunction, n, this};

//     gsl_vector* x = gsl_vector_alloc(n);
//     gsl_vector_set(x, 0, root_guess[0]); // mu_e
//     gsl_vector_set(x, 1, root_guess[1]); // mu_p

//     const gsl_multiroot_fsolver_type* T;
//     gsl_multiroot_fsolver* s;

//     T = gsl_multiroot_fsolver_hybrids;
//     s = gsl_multiroot_fsolver_alloc(T, n);
//     gsl_multiroot_fsolver_set(s, &f, x);

//     size_t iter = 0;
//     int status;
//     const double tol = 1e-8;

//     do {
//         iter++;
//         status = gsl_multiroot_fsolver_iterate(s);

//         if (status)
//             break;

//         status = gsl_multiroot_test_residual(s->f, tol);
//     } while (status == GSL_CONTINUE && iter < 1000);

//     // Update root_guess with the solution
//     root_guess[0] = gsl_vector_get(s->x, 0);
//     root_guess[1] = gsl_vector_get(s->x, 1);

//     gsl_multiroot_fsolver_free(s);
//     gsl_vector_free(x);

//     return status;
// }

//--------------------------------------------------------------
// Solve the equilibrium equations
int Fermi_Gas_Many::SolveEquations()
{
    const size_t n = 1;

    gsl_multiroot_function f = {&EquilibriumFunction, n, this};

    gsl_vector* x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, root_guess[0]); // rho_n
    // gsl_vector_set(x, 1, root_guess[1]); // mu_p

    const gsl_multiroot_fsolver_type* T;
    gsl_multiroot_fsolver* s;

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, n);
    gsl_multiroot_fsolver_set(s, &f, x);

    size_t iter = 0;
    int status;
    const double tol = 1e-8;

    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        if (status)
            break;

        status = gsl_multiroot_test_residual(s->f, tol);
    } while (status == GSL_CONTINUE && iter < 1000);

    // Update root_guess with the solution
    root_guess[0] = gsl_vector_get(s->x, 0);
    std::cout <<  root_guess[0] << "\n" ;
    // root_guess[1] = gsl_vector_get(s->x, 1);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return status;
}


//--------------------------------------------------------------
// Compute energy density
double Fermi_Gas_Many::EDens(const double& in_rho)
{
    // Set the total baryon density
    rho = in_rho;

    // Provide an initial guess for the solver
    FindGuess();

    // Solve the equilibrium equations
    int status = SolveEquations();

    if (status != GSL_SUCCESS)
    {
      std::cerr << "Solver failed to converge." << std::endl;
      return 0.0;
    }

    // Compute total energy density
    double total_energy_density = 0.0;

    for (const auto& baryon : baryon_pool)
    {
      total_energy_density += EDens_Individual(baryon);
    }

    for (const auto& lepton : lepton_pool)
    {
      total_energy_density += EDens_Individual(lepton);
    }

    return total_energy_density; // Units: MeV/fm^3 or appropriate
}


//--------------------------------------------------------------
// Energy density in g/cm^3 
// --> input rho has units of fm^{-3}
double Fermi_Gas_Many::EDens_Individual(Particle* in_p) 
{
  double rho_p = in_p->GetRho() ;

  // kF has units of fm^{-1}
  double kF = pow(3*M_PI*M_PI*rho_p, 1.0/3.0) ;

  // Converting the mass from MeV to fm^{-1}
  double m2 = in_p->M*MEV_2_INV_FM ;

  double mu = sqrt( m2*m2 + kF*kF ) ;

  // eps has unit of fm^{-4}
  double eps = mu*kF*(mu*mu - m2*m2/2) ;
        eps -= (pow(m2, 4)/2) * log((mu + kF) / m2) ;
        eps /= 4*M_PI*M_PI ;

  // converting eps unit from fm^{-4} to g/cm^3
  return eps*INV_FM4_2_G_CM3 ;
}

//--------------------------------------------------------------
// Compute pressure
double Fermi_Gas_Many::Press(const double& in_rho)
{
    // Ensure equilibrium equations are solved
    EDens(in_rho);

    // Compute total pressure
    double total_pressure = 0.0;

    for (const auto& baryon : baryon_pool)
    {
        total_pressure += baryon->PTerm();
    }

    for (const auto& lepton : lepton_pool)
    {
        total_pressure += lepton->PTerm();
    }

    return total_pressure; // Units: MeV/fm^3 or appropriate
}

//--------------------------------------------------------------
// Generate EOS data row
std::vector<double> Fermi_Gas_Many::EOSRow(const double& in_rho)
{
    std::vector<double> row;

    double energy_density = EDens(in_rho);
    double pressure = Press(in_rho);

    row.push_back(in_rho);
    row.push_back(energy_density);
    row.push_back(pressure);

    // Add particle fractions
    for (const auto& baryon : baryon_pool)
    {
        double fraction = baryon->GetRho() / in_rho;
        row.push_back(fraction);
    }

    for (const auto& lepton : lepton_pool)
    {
        double fraction = lepton->GetRho() / in_rho;
        row.push_back(fraction);
    }

    return row;
}

//--------------------------------------------------------------
// Provide EOS data header
std::string Fermi_Gas_Many::EOSHeader() const
{
    std::string header = "rho\tEnergyDensity\tPressure\t";

    for (const auto& baryon : baryon_pool)
    {
        header += baryon->GetName() + "_Frac\t";
    }

    for (const auto& lepton : lepton_pool)
    {
        header += lepton->GetName() + "_Frac\t";
    }

    return header;
}
//--------------------------------------------------------------
