#include "tests/relativity/candidate_capture.hpp"
// -*- lsst-c++ -*-
// Copyright (c) 2026 Mohammadreza Zakeri. MIT License — see LICENSE.
/**
 * Phase 4D-BL: production regression AFTER independent verification (ADR-0007/0008/0009).
 * This artifact is not an independent physics oracle. See PHASE4D_MONOPOLE_BASELINE.md.
 * Normal execution only compares; --emit requires a fresh path outside tests/baselines.
 */
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

namespace fs = std::filesystem;
using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;
using CompactStar::Core::TOVSolveStatus;

namespace
{
// Expose the existing cutoff for inspection only; no second surface locator or formula.
struct CheckedTOV : TOVSolver
{
    using TOVSolver::PressureCutoff;
};

using Row = std::array<double, 14>;
constexpr std::array<const char *, 14> kColumns = {
    "target_M_Msun", "epsilon_c_gcm3", "achieved_M_Msun", "R_km",
    "p_cut_dyn_cm2", "epsilon_surface_gcm3", "I_km3",
    "m0_surface_over_Omega2_km3", "p0star_surface_over_Omega2_km2",
    "delta_p0_surface_over_Omega2_dimensionless", "xi0_surface_over_Omega2_km3",
    "surface_shell_mass_over_Omega2_km3", "I2_over_R3_km3", "deltaM_over_Omega2_km3"};

void Require(bool ok, const std::string &message)
{
    if (!ok)
        throw std::runtime_error(message);
}

// Arithmetic identity bound already used by hartle_monopole_cmf; NOT a regression tolerance.
void Identity(double got, long double expected, const char *name)
{
    Require(std::isfinite(got) && std::isfinite(expected) &&
                std::abs(static_cast<long double>(got) - expected) <=
                    1e-14L * std::abs(expected), name);
}

struct WorkDirectory
{
    fs::path path;
    WorkDirectory()
    {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        for (int i = 0; i < 100; ++i)
        {
            path = fs::temp_directory_path() /
                ("compactstar_monopole_regression_" + std::to_string(stamp) + "_" + std::to_string(i));
            if (fs::create_directory(path))
                return;
        }
        throw std::runtime_error("cannot create private producer working directory");
    }
    ~WorkDirectory()
    {
        std::error_code ec;
        fs::remove_all(path, ec);
    }
};

void CheckEmitPath(const fs::path &path)
{
    const auto resolved = fs::weakly_canonical(path);
    const auto baseline_dir = fs::canonical(COMPACTSTAR_BASELINE_DIR);
    std::string previous;
    fs::path prefix;
    for (const auto &part : resolved)
    {
        prefix /= part;
        Require(prefix != baseline_dir && !(previous == "tests" && part == "baselines"),
                "--emit must be outside tests/baselines (including symlinks)");
        previous = part.string();
    }
    Require(!fs::exists(path), "--emit requires a fresh candidate path");
    Require(fs::is_directory(resolved.parent_path()), "candidate parent directory must exist");
}

std::vector<Row> Compute(const fs::path &cold)
{
    WorkDirectory work;
    std::vector<Row> rows;
    for (const double target : {1.0, 1.4, 1.6, 2.0})
    {
        CheckedTOV tov; // Inherited production/default radial_res = 10000; no overrides.
        tov.SetWrkDir(work.path.string());
        tov.ImportEOS(cold.string(), true);
        std::vector<TOVPoint> points;
        std::vector<std::string> labels;
        const int count = tov.SolveToProfile(target, points, &labels);
        Require(count > 0 && points.size() >= 4 &&
                    tov.LastSolveStatus() == TOVSolveStatus::SURFACE_REACHED,
                "target did not return a complete SURFACE_REACHED profile");
        const double cutoff = tov.PressureCutoff();
        Require(points.back().p == cutoff, "final TOV pressure differs from governed p_cut");
        Require(std::abs(points.back().m - target) < 1e-4, "production target-mass contract failed");

        NStar star(points, labels);
        const NStar &current = star; // Read-only access must not bump Profile::Version().
        const auto &profile = current.Profile();
        const auto &sequence = current.GetSequence();
        Require(current.MonopoleResponse() == nullptr, "construction implicitly computed monopole");
        const auto &first = current.RotationResponse();
        Require(first.valid && profile.HasEosDEdP(), "first-order response or EOS derivative absent");
        Require(star.ComputeHartleMonopoleResponse(), "explicit production monopole solve failed");
        const auto *response = current.MonopoleResponse();
        Require(response && response->MatchesSource(&profile, profile.Version()),
                "absent, invalid or stale response refused");
        const int last = static_cast<int>(profile.GetRadius()->Size()) - 1;
        const auto size = profile.GetRadius()->Size();
        Require(response->r_grid == profile.GetRadius() &&
                    response->m0_over_Omega2.Size() == size &&
                    response->p0star_over_Omega2.Size() == size &&
                    response->delta_p0_over_Omega2.Size() == size &&
                    response->xi0_over_Omega2.Size() == size,
                "response does not cover its source profile");
        Require(response->R_surface == points.back().r && response->R_surface == sequence.r &&
                    response->I == first.I && std::abs(sequence.m - target) < 1e-4,
                "response/background or target-mass mismatch");

        const double radius = response->R_surface;
        const double exterior = first.I * first.I / (radius * radius * radius);
        const double epsilon = (*profile.GetEnergyDensity())[last]; // geometric km^-2
        const double pressure = (*profile.GetPressure())[last];     // geometric km^-2
        Identity(response->surface_shell_mass_over_Omega2,
                 4.0L * std::acos(-1.0L) * radius * radius * epsilon *
                     response->surface_xi0_over_Omega2, "surface shell identity failed");
        Identity(response->delta_p0_over_Omega2[last],
                 (static_cast<long double>(epsilon) + pressure) * response->p0star_over_Omega2[last],
                 "Eulerian pressure coefficient identity failed");
        Require(response->surface_xi0_over_Omega2 == response->xi0_over_Omega2[last],
                "surface displacement differs from terminal coefficient");
        const long double sum = static_cast<long double>(response->m0_over_Omega2[last]) +
                                response->surface_shell_mass_over_Omega2 + exterior;
        Identity(response->deltaM_over_Omega2, sum, "deltaM decomposition failed");

        Row row = {target, sequence.ec, sequence.m, radius, cutoff, points.back().e, first.I,
                   response->m0_over_Omega2[last], response->p0star_over_Omega2[last],
                   response->delta_p0_over_Omega2[last], response->surface_xi0_over_Omega2,
                   response->surface_shell_mass_over_Omega2, exterior, response->deltaM_over_Omega2};
        for (std::size_t i = 0; i < row.size(); ++i)
            Require(std::isfinite(row[i]) && row[i] > 0.0,
                    std::string("nonfinite or wrong-sign canonical value: ") + kColumns[i]);
        // Reuse the public provenance guard at the point the response enters the artifact.
        Require(current.MonopoleResponse() == response &&
                    response->MatchesSource(&profile, profile.Version()), "response became stale");
        rows.push_back(row);
        std::cout << std::setprecision(17) << "target=" << target
                  << " SURFACE_REACHED p=p_cut current=true residual_M=" << std::abs(sequence.m - target)
                  << " decomposition_residual_km3=" << static_cast<long double>(row.back()) - sum << '\n';
    }
    return rows;
}

std::string Serialize(const std::vector<Row> &rows)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << "# CompactStar Hartle monopole regression; schema=1; Phase 4D-BL\n"
           "# Independently verified before this artifact; regression output, not a physics oracle.\n"
           "# Verification anchor: eccbfa6951ec7ed489e7dfde1fc93c2759d57e2a (4D-RV + 4D-DA).\n"
           "# EOS: DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos; production SolveToProfile.\n"
           "# Default radial_res=10000; targets 1.0/1.4/1.6/2.0 Msun; mass_tol=1e-4 Msun absolute.\n"
           "# R_* is the finite p=p_cut event, not the p=0 vacuum surface (ADR-0009).\n"
           "# All _over_Omega2 fields are per Omega_geom^2; Omega_geom=Omega_phys/c [km^-1].\n"
           "# Explicit production monopole solve at fixed selected epsilon_c; no implicit physical spin.\n"
           "# deltaM_hat = m0_hat(R_*^-) + surface_shell_hat + I^2/R_*^3 (ADR-0007/0008).\n"
           "# 17 significant digits, classic locale, round-trip doubles; deterministic byte comparison.\n";
    for (std::size_t i = 0; i < kColumns.size(); ++i)
        out << (i ? "\t" : "") << kColumns[i];
    out << '\n' << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto &row : rows)
    {
        for (std::size_t i = 0; i < row.size(); ++i)
            out << (i ? "\t" : "") << row[i];
        out << '\n';
    }
    return out.str();
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        Require(argc == 2 || (argc == 4 && (std::string(argv[2]) == "--emit" ||
                                            std::string(argv[2]) == "--compare")),
                "usage: hartle_monopole_regression <EOS_DATA_ROOT> [--emit <fresh-tsv> | --compare <tsv>]");
        const bool emit = argc == 4 && std::string(argv[2]) == "--emit";
        const fs::path path = argc == 4 ? fs::path(argv[3]) :
            fs::path(COMPACTSTAR_BASELINE_DIR) / "hartle_monopole_dscmf1_debug.tsv";
        std::string expected;
        if (emit)
            CheckEmitPath(path);
        else
        {
            std::ifstream in(path, std::ios::binary);
            Require(in.is_open(), "baseline missing: " + path.string());
            expected.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
            Require(!in.bad(), "baseline read failed");
        }
        const auto cold = fs::path(argv[1]) / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
        Require(fs::is_regular_file(cold), "required EOS missing: " + cold.string());
        const std::string bytes = Serialize(Compute(cold));
        unit_candidate_evidence::Capture("hartle_monopole_dscmf1_debug.tsv", [&](const fs::path &p) {
            std::ofstream out(p, std::ios::binary); out.exceptions(std::ios::failbit | std::ios::badbit); out << bytes; out.close();
        });
        if (emit)
        {
            std::ofstream out;
            out.exceptions(std::ios::failbit | std::ios::badbit);
            out.open(path, std::ios::binary);
            out << bytes;
            out.close();
            std::cout << "emitted 4 canonical rows to " << path << '\n';
        }
        else
        {
            Require(bytes == expected, "regression bytes differ from " + path.string());
            std::cout << "PASS: 4 canonical rows reproduce the regression artifact byte-for-byte\n";
        }
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
