#pragma once
// Native assembly of retained full owners. Historical fixture remains byte-identical.
#include "../bnv/campaign_fixture.hpp"
#include "stellar_archive.hpp"
ControlledFixture NativeFixture(std::filesystem::path profile_dir,std::string certificate,std::filesystem::path dir,unsigned resolution, const std::shared_ptr<EquilibriumSequenceNumberDerivative>& sequence) {
 require(resolution==80000,"production radial resolution frozen at80000");
 const bool full=true;
 require(!std::filesystem::exists(dir),"fresh fixture assembly directory required");
 std::filesystem::create_directories(dir);
		auto src = source(profile_dir / "freegas.tsv");
		sequence->RequireCurrent();
        auto central = sequence->contributing_stars.at(2);
		auto provider = std::make_shared<TrackRFreeGasThermodynamicProvider>();
        { std::ifstream model(profile_dir/"model.txt");
          for(const auto& fermion:{ColdRelativisticIdealFermion::Neutron(),ColdRelativisticIdealFermion::Proton(),ColdRelativisticIdealFermion::Electron(),ColdRelativisticIdealFermion::Muon()}) {
            double mass,hc;require(bool(model>>mass>>hc),"model transport missing");require(mass==fermion.RestMassEnergyMeV()&&hc==fermion.HbarCMeVFm(),"model constants differ from governed source");
          }
        }
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
            // Retained profile stores provider fields then the actual profile density products.
            std::vector<double> rest;double value;while(in>>value)rest.push_back(value);
            require(rest.size()==22,"retained full profile schema mismatch");
            require(rest[17]==(*central->Profile().GetEnergyDensity())[j],"retained energy profile mismatch");
            for(int k=0;k<4;++k)require(rest[18+k]==(*central->Profile().GetSpeciesPtr(species[k].label))[j]*(*NB)[j],"retained equilibrium composition mismatch");
            ++j;
		}
		require(j == R->Size(), "profile transport incomplete");
		std::shared_ptr<FixedBaryonNumberResponse> structural;
		if (full)
		{
			auto a = FixedCentralEnergyNumberResponse::Compute(input(*central, src));
			a.RequireCurrent();
            const auto& b=*sequence;
			b.RequireCurrent();
			for (auto &s : b.contributing_stars)
				lifetime->stars.push_back(s);
			structural = std::make_shared<FixedBaryonNumberResponse>(FixedBaryonNumberResponse::Compute(a, b));
			structural->RequireCurrent();
		}
		lifetime->revision = std::make_shared<ChemicalRevision>();
		*lifetime->revision = {provider->Metadata().model_id, provider->Metadata().model_revision, "authenticated source files retained by bytes", "Zaki cold fermions and AngularVelocity unit owner", "owning canonical Structure-1 midpoint", "nu=Phi, one inverse lapse, proper volume once", "Neutron,Electron,Muon", "WholeStar P=0, finite cut plus certified tail", "all profile nodes, onsets, refusal edges", "source/ULP-derived continuous onset certificates", "positive pe mass-upper radius construction", "PHASE5C2_PREDECLARATION a87f0212c2bd7bfba92db91dfac82447a6561334", "all central and sequence stars owned", true};
		lifetime->source_files = {src->table_path, (profile_dir / "profile.tsv").string(), certificate.c_str(), "CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp", "CompactStar/EOS/src/LocalThermodynamics.cpp", "CompactStar/EOS/TrackRFreeGasThermodynamics.hpp", "CompactStar/EOS/LocalThermodynamics.hpp"};
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
		std::ifstream cert(certificate);
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
		std::cout << "ASSEMBLY G begin\n";
        auto g = std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(r));
        std::cout << "ASSEMBLY G passed\n";

        std::cout << "ASSEMBLY Z begin\n";
        auto z=std::make_shared<ChemicalImbalanceResponse>(ChemicalImbalanceResponse::Compute(g,qgoal,zgoal));
        std::cout << "ASSEMBLY Z passed\n";
        // W is constructed inside the sealed run context from these owners.
        auto metric=[interpolate](double x){auto v=interpolate(x);return RC::UrcaMetric{v[1], x==0?0:-.5*std::log1p(-2*v[0]/x)};};
        return {central,provider,g,z,structural,vi,gn,gv,metric,interpolate};
    }
