// Diagnostic only: reuse the existing static solve helper, consume no certificate.
#define main phase5b_original_validation_main
#include "../../tests/analysis/phase5b_freegas_validation.cpp"
#undef main
int main(int argc, char** argv) {
    try {
        require(argc == 4, "profile-directory fresh-work output.tsv");
        const std::filesystem::path profile(argv[1]), work(argv[2]);
        require(!std::filesystem::exists(work), "fresh diagnostic directory required");
        std::filesystem::create_directories(work);
        auto src = source(profile / "freegas.tsv");
        auto star = solve(src->table_path, 1.10e15, 80000, work / "central");
        auto R = star->Profile().GetRadius(), M = star->Profile().GetMass();
        auto N = star->Profile().GetMetricNu(), B = star->Profile().GetBaryonDensity();
        std::ifstream original(profile / "profile.tsv");
        require(bool(original), "authenticated profile missing");
        std::ofstream output(argv[3]);
        require(bool(output), "cannot write diagnostic");
        output << std::setprecision(17) << "index\tr_mac\tr_linux\tm_mac\tm_linux\tnu_mac\tnu_linux\tnb_mac\tnb_linux\n";
        std::string line; std::getline(original, line);
        std::size_t j = 0;
        while (std::getline(original, line)) {
            std::istringstream in(line); double r, m, nu, nb;
            require(bool(in >> r >> m >> nu >> nb) && j < R->Size(), "diagnostic profile shape differs");
            output << j << '\t' << r << '\t' << (*R)[j] << '\t' << m << '\t' << (*M)[j]
                   << '\t' << nu << '\t' << (*N)[j] << '\t' << nb << '\t' << (*B)[j] << '\n';
            ++j;
        }
        require(j == R->Size(), "diagnostic profile row count differs");
        std::cout << "PROFILE_DIAGNOSTIC_ONLY rows " << j << " certificates_not_consumed\n";
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
