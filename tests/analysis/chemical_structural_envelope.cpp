// ADR-0013 / Phase-5C-1R M1 and M2, TEST SIDE ONLY.
// Reuse the authenticated Phase-5B fixture and its independent PB10 quadratures
// without editing a source whose bytes belong to the governed baseline producer.
// No chemical coefficient is computed by this executable.
#define main phase5b_original_validation_main
#include "phase5b_freegas_validation.cpp"
#undef main

int main(int argc, char **argv)
{
    try
    {
        require(argc == 2, "expected fresh evidence directory");
        const std::filesystem::path dir(argv[1]);
        require(!std::filesystem::exists(dir), "evidence directory must be fresh");
        std::filesystem::create_directories(dir);
        std::cout << std::setprecision(17) << std::unitbuf;
        std::ofstream raw(dir / "radial.tsv"), pb10(dir / "pb10.tsv");
        raw << std::setprecision(17)
            << "resolution\tspecies\tA\tB\tA_B\tB_B\tcentral_energy_per_q\tK\tK_numerical_error\tI_phys\tE_I_numerical\n";
        pb10 << std::setprecision(17)
             << "species\tK_prod\tK_PB10\tDelta_K\trelative_discrepancy\tunchanged_PB10_budget\n";
        const TrackRFreeGasThermodynamicProvider eos;
        const auto table = generate(eos, dir / "freegas.tsv", 8192);
        const auto src = source(table.file);
        // Protocol fixed in PHASE5C2_PREPRODUCTION_STRUCTURAL_ENVELOPE.md
        // before this executable's first run. Complete A/B/K at EVERY rung.
        for (std::size_t resolution : {20000, 40000, 80000})
        {
            const auto rung = dir / ("radial-" + std::to_string(resolution));
            auto central = solve(table.file.string(), 1.10e15, resolution, rung / "central");
            auto a = FixedCentralEnergyNumberResponse::Compute(input(*central, src));
            a.RequireCurrent();
            NumberSequenceRecipe recipe{src, species, whole, 1.10e15, 1.095e15, 1.105e15,
                "Structure-1 monotone midpoint smooth central branch", {.001, .0005, .00025},
                resolution, (rung / "sequence").string(), tail_bound};
            recipe.tail_policy_identity = "positive-source pe comparison inequalities";
            recipe.tail_policy_revision = "PB13-comparison-v1";
            auto b = EquilibriumSequenceNumberDerivative::Compute(recipe);
            b.RequireCurrent();
            auto k = FixedBaryonNumberResponse::Compute(a, b);
            k.RequireCurrent();
            require(b.contributing_stars.size() == 15, "complete 15-star sequence required");
            const auto av = a.Values(), bv = b.Values(), kv = k.Values(), ke = k.Errors();
            const auto physical = k.WholeStarIPhysical();
            const double inverse_c = AngularVelocity::FromRadPerSecond(1).GeomKmInverse();
            for (int i = 0; i < 4; ++i)
            {
                raw << resolution << '\t' << i << '\t' << av[i] << '\t' << bv[i]
                    << '\t' << k.A_B << '\t' << k.B_B << '\t' << k.central_energy_per_q
                    << '\t' << kv[i] << '\t' << ke[i] << '\t' << physical[i]
                    << '\t' << ke[i] * inverse_c * inverse_c << std::endl;
                std::cout << "M1 resolution=" << resolution << " species=" << i
                          << " A=" << av[i] << " B=" << bv[i] << " K=" << kv[i] << '\n';
            }
            if (resolution != 80000) continue;
            // Exact PB10 reconstruction and unchanged acceptance inequalities.
            // This is test mathematics, not another implementation of production K.
            const long double AB = static_cast<long double>(av[0]) + av[1];
            const long double BB = static_cast<long double>(bv[0]) + bv[1];
            const auto current0 = finite_current(*central, 0);
            const auto current1 = finite_current(*central, 1e-7);
            const auto current2 = finite_current(*central, 5e-8);
            std::array<double, 4> charged_K{};
            const double h = .00025, eps = b.metadata.central_energy_km_minus2;
            for (int i = 0; i < 4; ++i)
            {
                const double ac = 2 * (current2[i] - current0[i]) / 5e-8
                    - (current1[i] - current0[i]) / 1e-7;
                const double bc = (8 * (direct_nodal_count(*b.contributing_stars[13], i)
                    - direct_nodal_count(*b.contributing_stars[11], i))
                    - (direct_nodal_count(*b.contributing_stars[14], i)
                    - direct_nodal_count(*b.contributing_stars[10], i))) / (12 * h * eps);
                charged_K[i] = ac - bc * double(AB / BB);
                const double delta = charged_K[i] - kv[i];
                pb10 << i << '\t' << kv[i] << '\t' << charged_K[i] << '\t' << delta
                     << '\t' << std::abs(delta / kv[i]) << '\t' << ke[i] << std::endl;
                require(std::abs(delta) <= ke[i],
                    "PB10 independent charged-current/sequence comparison exceeds propagated budget");
            }
            require(std::abs(charged_K[1] - (charged_K[2] + charged_K[3])) <= k.charge_budget,
                    "PB10 independent neutral reconstruction failed");
            require(std::abs(kv[1] - kv[0] - kv[3]) > k.charge_budget,
                    "wrong charge map/order mutation escaped");
            std::cout << "M2 raw per-species PB10 reporting PASS\n";
        }
        require(bool(raw) && bool(pb10), "raw evidence write failed");
        std::cout << "M1 measurement completed; envelope adjudication remains separate\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "STOP " << error.what() << '\n';
        return 1;
    }
}
