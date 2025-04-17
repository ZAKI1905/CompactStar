/*
  Fermi_Gas_Many Model Class Header
  Author: [Your Name]
  Date: [Current Date]
  Description:
    This class models the equation of state (EOS) for neutron star matter
    by including multiple fermion species (protons, neutrons, electrons, muons).
    It calculates the energy density and pressure by finding the number densities
    of each particle species that minimize the total energy under the constraints
    of electric charge neutrality and baryon number conservation.
*/

#ifndef CompactStar_Fermi_Gas_Many_HPP
#define CompactStar_Fermi_Gas_Many_HPP

#include <vector>
#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "Model.hpp"
#include "Particle.hpp" // Base class for Lepton and Baryon
#include "Baryon.hpp"
#include "Lepton.hpp"

namespace CompactStar {

class Fermi_Gas_Many : public Model {
private:
    // Particle instances
    Baryon proton;
    Baryon neutron;
    Lepton electron;
    // Lepton muon;

    // Particle pools
    std::vector<Particle*> baryon_pool;
    std::vector<Particle*> lepton_pool;

    // Total baryon density
    double rho = 0 ;

    // Solver variables
    size_t eq_size;
    std::vector<double> root_guess;
    std::vector<double> equilibrium_eqs;

    // Methods
    void SetEquilibriumEqs();
    double chem_equilibrium_eq ;
    void UpdateParticles();
    void FindGuess();

    // GSL solver functions
    static int EquilibriumFunction(const gsl_vector* x, void* params, gsl_vector* f);

    // Solver
    int SolveEquations();

public:
    // Constructor and Destructor
    Fermi_Gas_Many();
    ~Fermi_Gas_Many();

    // EOS computation
    double EDens(const double& in_rho) override;
    double Press(const double& in_rho) override;

    // Energy density of individual particles in g/cm^3 
    // --> input rho has units of fm^{-3}
    double EDens_Individual(Particle* in_p) ;

    // Data output
    std::vector<double> EOSRow(const double& in_rho) override;
    std::string EOSHeader() const override;
};

} // namespace CompactStar

#endif // CompactStar_Fermi_Gas_Many_HPP