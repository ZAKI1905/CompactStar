// TEST-ONLY fresh canonical Track-R producer. No chemical integration here.
#include <CompactStar/EOS/TrackRFreeGasThermodynamics.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/Core/NStar.hpp>
#include <CompactStar/RelativityUnits.hpp>
#include <CompactStar/Units.hpp>
#include <tests/eos/structure1/table.hpp>
#include "chemical_enthalpy_reference.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <type_traits>
#include <limits>
using namespace CompactStar;
namespace {
struct Solver : Core::TOVSolver {
    const auto& Labels() const { return eos_tab.extra_labels; }
};
int Dimension(const TrackRFreeGasThermodynamicProvider& eos, double n) {
    try {
        return std::visit([](const auto& v) {
            return int(std::decay_t<decltype(v)>::response_dimension);
        }, eos.EquilibriumAt(n));
    } catch (const EquilibriumResolutionError&) { return -1; }
}
}
int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("directory table-resolution radial-resolution");
        const std::filesystem::path dir(argv[1]);
        if (std::filesystem::exists(dir)) throw std::runtime_error("fresh output directory required");
        std::filesystem::create_directories(dir);
        TrackRFreeGasThermodynamicProvider eos;
        const auto table = structure1::generate(eos, dir/"freegas.tsv", std::stoul(argv[2]));
        Solver solver;
        solver.SetWrkDir((dir/"star").string());
        solver.ImportEOS(table.file.string(), true);
        solver.SetRadialRes(std::stoul(argv[3]));
        std::vector<Core::TOVPoint> points;
        if (solver.SingleStarSolveToTOVPoints(1.10e15, points) <= 0 ||
            solver.LastSolveStatus() != Core::TOVSolveStatus::SURFACE_REACHED)
            throw std::runtime_error("canonical TOV failed to reach finite-cut surface");
        Core::NStar star(points, solver.Labels());
        const auto& profile = star.Profile();
        auto r=profile.GetRadius(), m=profile.GetMass(), nu=profile.GetMetricNu();
        auto nb=profile.GetBaryonDensity(), pressure=profile.GetPressure();
        std::ofstream out(dir/"profile.tsv");
        out << std::setprecision(17)
            << "r\tm\tnu\tnB\tPgeom\tdim\tnn\tnp\tne\tnmu\tepsMeVfm3\tPMeVfm3"
            << "\tH00\tH01\tH02\tH10\tH11\tH12\tH20\tH21\tH22"
            << "\teps_profile_geom\tnn_profile\tnp_profile\tne_profile\tnmu_profile\n";
        for (size_t i=0; i<r->Size(); ++i) {
            // Padding is transport serialization ONLY; Python slices the active H.
            std::array<std::array<double,3>,3> h{};
            int dimension=-1;
            const auto value=eos.BarotropeAt((*nb)[i]);
            try {
                std::visit([&](const auto& v) {
                    using V=std::decay_t<decltype(v)>;
                    dimension=V::response_dimension;
                    if constexpr (V::response_dimension>0)
                        for (int j=0;j<dimension;++j) for (int k=0;k<dimension;++k)
                            h[j][k]=v.hessian(j,k);
                }, eos.EquilibriumAt((*nb)[i]));
            } catch (const EquilibriumResolutionError&) {}
            out << (*r)[i] << '\t' << (*m)[i] << '\t' << (*nu)[i] << '\t'
                << (*nb)[i] << '\t' << (*pressure)[i] << '\t' << dimension;
            for (double n:value.number_densities_fm3) out << '\t' << n;
            out << '\t' << value.energy_density_MeV_fm3 << '\t' << value.pressure_MeV_fm3;
            for (auto row:h) for (double entry:row) out << '\t' << entry;
            out << '\t' << (*profile.GetEnergyDensity())[i];
            for (const char* label:{"10","11","0","1"}) {
                const auto* column=profile.GetSpeciesPtr(label);
                if (!column) throw std::runtime_error("missing authenticated species label");
                out << '\t' << (*column)[i]*(*nb)[i];
            }
            out << '\n';
        }
        std::ofstream meta(dir/"model.txt");
        meta << std::setprecision(17);
        for (auto f:{ColdRelativisticIdealFermion::Neutron(),ColdRelativisticIdealFermion::Proton(),
                     ColdRelativisticIdealFermion::Electron(),ColdRelativisticIdealFermion::Muon()})
            meta << f.RestMassEnergyMeV() << ' ' << f.HbarCMeVFm() << '\n';
        meta << eos.NeutronOnsetBaryonDensityFm3() << ' ' << eos.MuonOnsetBaryonDensityFm3() << '\n';
        meta << RelativityUnits::PressureDynCm2ToKmMinus2(Units::MEV_FM3_TO_ERG_CM3) << '\n';
        meta << eos.Metadata().model_id << '\n' << eos.Metadata().model_revision << '\n';
        // Empirical boundary brackets, including BOTH sides of muon onset.
        // The Python model bound encloses entire intervals, not just sample nodes.
        std::ofstream windows(dir/"windows.tsv");
        windows << std::setprecision(17) << "onset\tside\tlast_refused\tfirst_available\n";
        for (int index=0;index<2;++index) {
            const double onset=index?eos.MuonOnsetBaryonDensityFm3():eos.NeutronOnsetBaryonDensityFm3();
            if (Dimension(eos,onset)!=0) throw std::runtime_error("threshold is not value-only");
            for (int side:{-1,1}) {
                if (index==0 && side<0) continue;
                double inner=0, outer=onset*(index?1e-8:1e-4);
                if (Dimension(eos,onset+side*outer)<=0) throw std::runtime_error("no available exterior bracket");
                for (int it=0;it<100;++it) {
                    double mid=inner+(outer-inner)/2;
                    if (mid==inner || mid==outer) break;
                    if (Dimension(eos,onset+side*mid)>0) outer=mid; else inner=mid;
                }
                windows << onset << '\t' << side << '\t' << onset+side*inner << '\t' << onset+side*outer << '\n';
            }
        }
        if (!out || !meta || !windows) throw std::runtime_error("fixture output failed");
        if (std::stoul(argv[2])==8192 && std::stoul(argv[3])==80000) {
            chemical_reference::Reference reference;
            std::ofstream ref(dir/"independent-star.tsv");
            ref << std::setprecision(17) << "R\tM\tGnn\tGee\tGem\tGmm\tBn\tBp\tBe\tBmu\n";
            for (const auto controls:{std::array<double,2>{1e-10,1e-8},
                                       std::array<double,2>{1e-12,1e-10}}) {
                const auto answer=reference.compute(controls[0],controls[1]);
                for (size_t j=0;j<answer.size();++j) ref << (j?"\t":"") << answer[j];
                ref << '\n';
            }
            if (!ref) throw std::runtime_error("independent reference output failed");
        }
        std::cout << "fresh profile nodes=" << r->Size() << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
