// TEST-ONLY governed fixture assembly. Numerical chemical work is production.
#define main phase5b_original_validation_main
#include "phase5b_freegas_validation.cpp"
#undef main
#include <CompactStar/Analysis/ChemicalResponse.hpp>
#include <gsl/gsl_version.h>
using M = ChemicalMatrix;
M matrix(std::istream &in, int n)
{
	M a(n, std::vector<double>(n));
	for (auto &r : a)
		for (auto &x : r)
			require(bool(in >> x), "incomplete certificate matrix");
	return a;
}
void emit_matrix(const char *label, const M &a)
{
	std::cout << "RESULT " << label << ' ' << a.size();
	for (const auto &r : a)
		for (double v : r)
			std::cout << ' ' << v;
	std::cout << '\n';
}
void emit_vector(const char *label, const std::vector<double> &a)
{
	std::cout << "RESULT " << label << ' ' << a.size();
	for (double v : a)
		std::cout << ' ' << v;
	std::cout << '\n';
}
int main(int argc, char **argv)
{
	try
	{
		require(argc == 6, "profile-directory certificate-path output-directory radial-resolution full-or-global");
		std::cout << std::setprecision(17) << std::unitbuf;
		std::cout << "PROVENANCE compiler " << __VERSION__ << "\nPROVENANCE gsl " << gsl_version << "\nPROVENANCE cplusplus " << __cplusplus << '\n';
#ifdef NDEBUG
		std::cout << "PROVENANCE configuration NDEBUG\n";
#else
		std::cout << "PROVENANCE configuration assertions-enabled\n";
#endif

		std::filesystem::path profile_dir(argv[1]), dir(argv[3]);
		require(!std::filesystem::exists(dir), "fresh production fixture directory required");
		std::filesystem::create_directories(dir);
		const auto resolution = std::stoul(argv[4]);
		const bool full = std::string(argv[5]) == "full";
		auto src = source(profile_dir / "freegas.tsv");
		auto central = solve(src->table_path, 1.10e15, resolution, dir / "central");
		auto provider = std::make_shared<TrackRFreeGasThermodynamicProvider>();
		auto lifetime = std::make_shared<ChemicalLifetime>();
		lifetime->provider = provider;
		lifetime->stars.push_back(central);
		// The independently characterized fresh profile must equal this owning star,
		// entry for entry, before its certificates are eligible for consumption.
		auto R = central->Profile().GetRadius(), Mass = central->Profile().GetMass(), Nu = central->Profile().GetMetricNu(), NB = central->Profile().GetBaryonDensity();
		std::ifstream profile(profile_dir / "profile.tsv");
		require(bool(profile), "fresh characterized profile missing");
		std::string line;
		std::getline(profile, line);
		std::size_t j = 0;
		while (std::getline(profile, line))
		{
			std::istringstream in(line);
			double r, m, nu, nb;
			require(bool(in >> r >> m >> nu >> nb) && j < R->Size(), "profile transport shape");
			require(r == (*R)[j] && m == (*Mass)[j] && nu == (*Nu)[j] && nb == (*NB)[j], "characterization profile differs from owning canonical star");
			++j;
		}
		require(j == R->Size(), "profile transport incomplete");
		std::shared_ptr<FixedBaryonNumberResponse> structural;
		if (full)
		{
			auto a = FixedCentralEnergyNumberResponse::Compute(input(*central, src));
			a.RequireCurrent();
			NumberSequenceRecipe recipe{src, species, whole, 1.10e15, 1.095e15, 1.105e15, "Structure-1 monotone midpoint smooth central branch", {.001, .0005, .00025}, resolution, (dir / "sequence").string(), tail_bound};
			recipe.tail_policy_identity = "positive-source pe comparison inequalities";
			recipe.tail_policy_revision = "PB13-comparison-v1";
			auto b = EquilibriumSequenceNumberDerivative::Compute(recipe);
			b.RequireCurrent();
			for (auto &s : b.contributing_stars)
				lifetime->stars.push_back(s);
			structural = std::make_shared<FixedBaryonNumberResponse>(FixedBaryonNumberResponse::Compute(a, b));
			structural->RequireCurrent();
		}
		lifetime->revision = std::make_shared<ChemicalRevision>();
		*lifetime->revision = {provider->Metadata().model_id, provider->Metadata().model_revision, "authenticated source files retained by bytes", "Zaki cold fermions and AngularVelocity unit owner", "owning canonical Structure-1 midpoint", "nu=Phi, one inverse lapse, proper volume once", "Neutron,Electron,Muon", "WholeStar P=0, finite cut plus certified tail", "all profile nodes, onsets, refusal edges", "source/ULP-derived continuous onset certificates", "positive pe mass-upper radius construction", "PHASE5C2_PREDECLARATION a87f0212c2bd7bfba92db91dfac82447a6561334", "all central and sequence stars owned", true};
		lifetime->source_files = {src->table_path, (profile_dir / "profile.tsv").string(), argv[2], "CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp", "CompactStar/EOS/src/LocalThermodynamics.cpp", "CompactStar/EOS/TrackRFreeGasThermodynamics.hpp", "CompactStar/EOS/LocalThermodynamics.hpp"};
		const auto provider_copy = (dir / "authenticated-provider-source.cpp").string();
		std::filesystem::copy_file("CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp", provider_copy);
		lifetime->source_files.push_back(provider_copy);
		lifetime->source_files.push_back((profile_dir / "model.txt").string());
		auto interpolate = [central, R, Mass, Nu, NB](double r)
		{std::size_t k=std::upper_bound(R->Values().begin(),R->Values().end(),r)-R->Values().begin();
  require(r>=0&&r<=(*R)[-1],"chemical node outside owning profile");if(k==0){double t=r/(*R)[0];return std::array<double,3>{(*Mass)[0]*t*t*t,(*Nu)[0],(*NB)[0]};}
  if(k==R->Size())return std::array<double,3>{(*Mass)[-1],(*Nu)[-1],(*NB)[-1]};double t=(r-(*R)[k-1])/((*R)[k]-(*R)[k-1]);return std::array<double,3>{(*Mass)[k-1]+t*((*Mass)[k]-(*Mass)[k-1]),(*Nu)[k-1]+t*((*Nu)[k]-(*Nu)[k-1]),(*NB)[k-1]+t*((*NB)[k]-(*NB)[k-1])}; };
		ChemicalIntegrationRequest r;
		r.lifetime = lifetime;
		r.background = [interpolate](double x)
		{auto v=interpolate(x);return ChemicalMetric{v[0],v[1]}; };
		r.local = [provider, interpolate](double x)
		{const auto v=provider->EquilibriumAt(interpolate(x)[2]);std::size_t n=std::visit([](const auto &a){return a.response_dimension;},v);require(n>0,"threshold H query prevented");
  // Provider/background characterization is supplied separately. Here H is
  // the declared binary matrix; actual factorization/congruence arithmetic
  // and its residual are owned and propagated by the production adapter.
  auto c=ChargeNeutralNumberSusceptibility::Compute(v,M(n,std::vector<double>(n)));return ChemicalLocalResponse{c.Values(),c.NumericalError(),c.Support()}; };
		std::ifstream cert(argv[2]);
		require(bool(cert), "certificate transport missing");
		std::size_t n;
		require(bool(cert >> n), "partition count missing");
		for (std::size_t i = 0; i < n; ++i)
		{
			ChemicalInterval s;
			require(bool(cert >> s.left >> s.right >> s.left_continuous_onset >> s.right_continuous_onset), "partition incomplete");
			r.intervals.push_back(s);
		}
		require(r.intervals.front().left == 0 && r.intervals.back().right == (*R)[-1], "whole profile support mismatch");
		for (double node : R->Values())
			require(std::any_of(r.intervals.begin(), r.intervals.end(), [&](const auto &s)
								{ return s.right == node; }),
					"original profile node missing from sealed partition");
		require(bool(cert >> n), "refusal count missing");
		for (std::size_t i = 0; i < n; ++i)
		{
			ChemicalRefusalCertificate f;
			require(bool(cert >> f.left >> f.right >> f.containing_left >> f.containing_right >> f.first_cell >> f.last_cell >> f.radius_upper >> f.mass_upper >> f.nu_lower), "refusal incomplete");
			f.branch = i == 0 ? "radial left npe; radial right pe" : "radial left npemu; radial right npe";
			f.availability_authority = "Track-R local-r3 source guards and monotone free-gas model; see certificate evidence";
			f.geometry_extrema_authority = "entire containing profile cells, linear metric, outward radius edges";
			f.numerical_error = matrix(cert, 3);
			r.refusals.push_back(f);
		}
		ChemicalPeTailCertificate tc;
		require(bool(cert >> tc.cut_radius >> tc.cut_mass >> tc.cut_nu >> tc.h_upper >> tc.epsilon_upper >> tc.bootstrap_radius_upper >> tc.total_mass_upper >> tc.radius_upper >> tc.susceptibility_upper >> tc.lapse_error_upper >> tc.response_upper), "pe tail certificate incomplete");
		tc.source_hypotheses = "authenticated Track-R pe positive monotone energy/pressure, whole-star cut";
		r.pe_tail = tc;
		r.center_error = matrix(cert, 3);
		r.tail_error = matrix(cert, 3);
		r.background_error = matrix(cert, 3);
		r.absolute_goal = matrix(cert, 3);
		auto qgoal = matrix(cert, 2), zgoal = matrix(cert, 2);
		std::vector<double> vi(2), gn(2), gv(2);
		for (auto *v : {&vi, &gn, &gv})
			for (auto &x : *v)
				require(bool(cert >> x), "dual W goals incomplete");
		r.center_authority = "regular center integrated, leading missing-term characterization";
		r.tail_authority = "source pe comparison with mass upper enclosure and same-cut independent shell checks";
		r.background_error_authority = "fresh table/radial independent background, local provider and equilibrium anchor characterization";
		r.structural_zeros = {{NumberAxis::Neutron, NumberAxis::Electron}, {NumberAxis::Electron, NumberAxis::Neutron}, {NumberAxis::Neutron, NumberAxis::Muon}, {NumberAxis::Muon, NumberAxis::Neutron}};
		r.structural_zero_authority = "independent free-gas neutral y-chart polynomial: neutron/lepton blocks decouple";
		auto g = std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(r));
		emit_matrix("G", g->Values());
		emit_matrix("E_G", g->NumericalError());
		emit_matrix("E_quadrature", g->QuadratureError());
		emit_vector("eigenvalues", g->Eigenvalues());
		std::cout << "RESULT G_condition 1 " << g->Condition() << "\nRESULT node_count 1 " << g->NodeCount() << '\n';
		const auto &diag = g->Diagnostics();
		std::vector<double> zeros, versions;
		for (const auto &entry : diag.structural_zeros)
		{
			zeros.push_back(int(entry.first));
			zeros.push_back(int(entry.second));
		}
		emit_vector("structural_zero_indices", zeros);
		for (const auto &star : lifetime->stars)
			versions.push_back(star->Profile().Version());
		emit_vector("profile_versions", versions);

		for (std::size_t i = 0; i < diag.validation_integrals.size(); ++i)
			emit_matrix(("G_ladder_" + std::to_string(i)).c_str(), diag.validation_integrals[i]);
		emit_matrix("E_local_arithmetic", diag.local_arithmetic_error);
		emit_matrix("E_center", diag.center_error);
		emit_matrix("E_tail", diag.tail_error);
		emit_matrix("E_background", diag.background_error);
		emit_matrix("E_refusal", diag.refusal_error);
		emit_vector("partition", g->Partition());
		std::vector<double> support;
		for (auto a : g->Support())
			support.push_back(int(a));
		emit_vector("support", support);
		std::cout << "RESULT validation_node_count 1 " << diag.validation_node_count << '\n';
		auto z = std::make_shared<ChemicalImbalanceResponse>(ChemicalImbalanceResponse::Compute(g, qgoal, zgoal));
		auto q = z->ReductionDiagnostic();
		emit_matrix("Q", q.q);
		emit_matrix("E_Q", q.numerical_error);
		emit_matrix("E_Schur_arithmetic", q.schur_arithmetic);
		emit_matrix("Z", z->Values());
		emit_matrix("E_Z", z->NumericalError());
		emit_matrix("E_global_solve", z->GlobalSolveError());
		emit_matrix("E_Z_arithmetic", z->ArithmeticError());
		std::cout << "RESULT rho 1 " << z->Rho() << '\n';
		if (full)
		{
			auto w = RotochemicalSpinDrive::Compute(z, structural, vi, gn, gv);
			emit_vector("I", w.IPhysical());
			emit_vector("E_I_numerical", w.INumericalError());
			emit_vector("V_I_validation", w.IValidationEnvelope());
			emit_vector("W", w.Values());
			emit_vector("E_W_numerical", w.NumericalError());
			emit_vector("V_W_validation", w.ValidationEnvelope());
			emit_vector("E_W_arithmetic", w.ArithmeticError());
			emit_vector("spin_action", w.Evaluate(2, -1));
			auto expect = [&](auto function)
			{bool failed=false;try{function();}catch(const std::runtime_error &){failed=true;}require(failed,"GC14 mutation returned stale numbers"); };
			auto saved = *lifetime->revision;
			for (auto member : {&ChemicalRevision::model, &ChemicalRevision::provider_revision, &ChemicalRevision::provider_bytes, &ChemicalRevision::particle_constants, &ChemicalRevision::background, &ChemicalRevision::metric, &ChemicalRevision::basis, &ChemicalRevision::domain, &ChemicalRevision::partition, &ChemicalRevision::onset, &ChemicalRevision::tail, &ChemicalRevision::accuracy, &ChemicalRevision::lifetime_token})
			{
				lifetime->revision.get()->*member += "mutation";
				expect([&]
					   { w.Values(); });
				*lifetime->revision = saved;
			}
			double old = structural->A_B;
			structural->A_B += std::abs(old) * 1e-6;
			expect([&]
				   { w.Values(); });
			structural->A_B = old;
			auto saved_source = *src;
			src->revision += "mutation";
			expect([&]
				   { w.Values(); });
			*src = saved_source;
			// Independent mutations of the actual central and sequence profile versions.
			for (std::size_t index : {std::size_t(0), std::size_t(1)})
			{
				auto &profile = const_cast<Core::StarProfile &>(lifetime->stars[index]->Profile());
				auto saved_profile = profile;
				profile.Touch();
				expect([&]
					   { w.Values(); });
				profile = saved_profile;
			}
			auto original_owner = lifetime->stars[0];
			lifetime->stars[0] = lifetime->stars[1];
			expect([&]
				   { w.Values(); });
			lifetime->stars[0] = original_owner;
			auto original_central = structural->metadata.sources[0].star;
			structural->metadata.sources[0].star = nullptr;
			expect([&]
				   { w.Values(); });
			structural->metadata.sources[0].star = original_central;
			auto original_sequence = structural->metadata.sources[1].star;
			structural->metadata.sources[1].star = nullptr;
			expect([&]
				   { w.Values(); });
			structural->metadata.sources[1].star = original_sequence;
			// Both substituted sources remain alive and individually current.
			// Ordered dependency identity, rather than mere owner coverage, must refuse.
			std::swap(structural->metadata.sources[0], structural->metadata.sources[1]);
			expect([&]
				   { w.Values(); });
			std::swap(structural->metadata.sources[0], structural->metadata.sources[1]);
			std::swap(structural->metadata.sources[1], structural->metadata.sources[2]);
			expect([&]
				   { w.Values(); });
			std::swap(structural->metadata.sources[1], structural->metadata.sources[2]);
			auto metadata_copy = structural->metadata;
			structural->metadata.species[2].charge = 0;
			expect([&]
				   { RotochemicalSpinDrive::Compute(z, structural, vi, gn, gv); });
			structural->metadata = metadata_copy;
			structural->metadata.q_normalization = "q=Omega_phys^2";
			expect([&]
				   { RotochemicalSpinDrive::Compute(z, structural, vi, gn, gv); });
			structural->metadata = metadata_copy;
			structural->metadata.profile_units = "foreign units";
			expect([&]
				   { RotochemicalSpinDrive::Compute(z, structural, vi, gn, gv); });
			structural->metadata = metadata_copy;
			structural->metadata.profile_node_counts[0] += 1;
			expect([&]
				   { w.Values(); });
			structural->metadata = metadata_copy;
			// A genuinely destroyed foreign source is rejected using ownership
			// identity before either its raw star or profile can be dereferenced.
			auto temporary = std::make_shared<Core::NStar>();
			structural->metadata.sources[1].star = temporary.get();
			temporary.reset();
			expect([&]
				   { w.Values(); });
			structural->metadata.sources[1].star = original_sequence;
			// Same label/path, different bytes: generated transport, provider copy and EOS files only.
			for (const auto &filename : {std::string(argv[2]), src->table_path, provider_copy})
			{
				std::ifstream input_file(filename, std::ios::binary);
				std::string contents(std::istreambuf_iterator<char>(input_file), {});
				input_file.close();
				{
					std::ofstream f(filename, std::ios::app);
					f << "\n# deliberate GC14 mutation\n";
				}
				expect([&]
					   { w.Values(); });
				{
					std::ofstream f(filename, std::ios::binary);
					f << contents;
				}
				(void)w.Values();
			}
			// Dual acceptance is independently executable, including a zero-goal mutant.
			expect([&]
				   { RotochemicalSpinDrive::Compute(z, structural, vi, {0, 0}, gv); });
			expect([&]
				   { RotochemicalSpinDrive::Compute(z, structural, vi, gn, {0, 0}); });
			auto flipped = w.Evaluate(2, 1), normal = w.Evaluate(2, -1);
			for (int i = 0; i < 2; ++i)
				require(flipped[i] == -normal[i], "spin sign mutant not separated");
			// Release caller handles while the result retains all sources; then remove
			// an owner from its explicit bundle and refuse before any raw dereference.
			central.reset();
			provider.reset();
			g.reset();
			z.reset();
			structural.reset();
			(void)w.Values();
			auto owners = lifetime->stars;
			lifetime->stars.clear();
			expect([&]
				   { w.Values(); });
			lifetime->stars = owners;
			(void)w.Values();
			std::cout << "RESULT lifetime_mutations 1 1\n";
		}
		std::cout << "PRODUCTION FIXTURE COMPLETE\n";
		return 0;
	}
	catch (const std::exception &e)
	{
		std::cerr << "STOP " << e.what() << '\n';
		return 1;
	}
}
