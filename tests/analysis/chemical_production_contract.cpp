#include <CompactStar/Analysis/ChemicalResponse.hpp>
#include <CompactStar/EOS/TrackRFreeGasThermodynamics.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
using namespace CompactStar;
using namespace CompactStar::Analysis;
using M = ChemicalMatrix;
static_assert(!std::is_copy_assignable_v<GlobalChemicalNumberResponse>);
static_assert(!std::is_move_assignable_v<GlobalChemicalNumberResponse>);
static_assert(!std::is_copy_assignable_v<ChemicalImbalanceResponse>);
static_assert(!std::is_move_assignable_v<ChemicalImbalanceResponse>);
static_assert(!std::is_invocable_v<decltype(&ChargeNeutralNumberSusceptibility::Compute), ChemicalMatrix, ChemicalMatrix>);
static_assert(!std::is_invocable_v<decltype(&ChargeNeutralNumberSusceptibility::Compute), ChargeNeutralNumberSusceptibility, ChemicalMatrix>);
M zeros(int n = 3) { return M(n, std::vector<double>(n)); }
void require(bool b, const char *s)
{
	if (!b)
		throw std::runtime_error(s);
}
template <class F>
void refuse(F f, const std::string &reason)
{
	try
	{
		f();
	}
	catch (const std::runtime_error &e)
	{
		require(std::string(e.what()).find(reason) != std::string::npos, e.what());
		return;
	}
	throw std::runtime_error("required refusal missing: " + reason);
}
std::shared_ptr<ChemicalLifetime> owner()
{
	auto o = std::make_shared<ChemicalLifetime>();
	o->revision = std::make_shared<ChemicalRevision>();
	auto &r = *o->revision;
	r = {"analytic fixture", "revision-1", "independent analytic source", "unit constants", "synthetic background", "nu=Phi, one inverse lapse", "n,e,mu", "whole analytic domain", "sealed partition", "declared continuous", "exact empty tail", "analytic absolute goal", "owned analytic token", true};
	return o;
}
ChemicalIntegrationRequest request(double radius = 1)
{
	ChemicalIntegrationRequest r;
	r.lifetime = owner();
	r.background = [](double)
	{ return ChemicalMetric{0, 0}; };
	r.local = [](double)
	{M c=zeros();for(int i=0;i<3;++i)c[i][i]=1e-54;return ChemicalLocalResponse{c,zeros(),{NumberAxis::Neutron,NumberAxis::Electron,NumberAxis::Muon}}; };
	r.intervals = {{0, radius}};
	r.center_error = r.tail_error = r.background_error = zeros();
	r.absolute_goal = M(3, std::vector<double>(3, 1e-8));
	r.center_authority = "regular analytic center integrated directly";
	r.tail_authority = "exact finite support";
	r.background_error_authority = "exact manufactured metric";
	return r;
}
ActiveLocalThermodynamicEvaluation full(const M &h)
{
	ChargeNeutralChemicalHessian H;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			H.value_MeV_fm3[i][j] = h[i][j];
	return LocalThermodynamicEvaluation{MakeChargeNeutralCompositionState({.2, .04, .02}), 0, {}, H};
}
void emit(const char *name, const M &a)
{
	std::cout << name;
	for (auto &row : a)
		for (double x : row)
			std::cout << '\t' << x;
	std::cout << '\n';
}
int main(int argc, char **argv)
{
	try
	{
		std::cout << std::setprecision(17);
		if (argc > 1 && std::string(argv[1]) == "gc9")
		{
			for (int variable : {0, 1})
				for (unsigned order : {16, 32, 64})
					for (int mutant = 0; mutant < 7; ++mutant)
					{
						auto r = request(12);
						r.order = order;
						r.absolute_goal = M(3, std::vector<double>(3, 1e-7));
						// Fixed geometric bisection ladder: resolve the coarse GL8 error estimator
						// without changing the independent analytic acceptance goal.
						r.intervals.clear();
						for (int cell = 0; cell < 8; ++cell)
							r.intervals.push_back({1.5 * cell, 1.5 * (cell + 1)});
						r.background = [mutant](double x)
						{double m=1.8*std::pow(x/12,3),nu=std::log((3*std::sqrt(.7)-std::sqrt(1-2*m/x))/2);
    if(mutant==1)nu=0;if(mutant==2)nu*=2;if(mutant==3)nu*=.5;if(mutant==4)m=0;if(mutant==5)nu=-nu;return ChemicalMetric{m,nu}; };
						r.local = [variable, mutant](double x)
						{M c=zeros();double f=(1+variable*x*x/144)*1e-54;if(mutant==6)f*=1-3.6*x*x/(12*12*12);for(int i=0;i<3;++i)c[i][i]=f;return ChemicalLocalResponse{c,zeros(),{NumberAxis::Neutron,NumberAxis::Electron,NumberAxis::Muon}}; };
						auto g = GlobalChemicalNumberResponse::Compute(r);
						std::cout << "GC9\t" << variable << '\t' << order << '\t' << mutant << '\t' << g.Values()[0][0] << '\n';
					}
			return 0;
		}
		M h{{4, -2, -3}, {-2, 11, 8}, {-3, 8, 16}};
		auto local = ChargeNeutralNumberSusceptibility::Compute(full(h), zeros());
		emit("GC3_6_C", local.Values());
		emit("GC3_6_E", local.NumericalError());
		// Positive high condition is permitted when the response remains resolved.
		auto hi = ChargeNeutralNumberSusceptibility::Compute(full({{1e-24, 0, 0}, {0, 1, 0}, {0, 0, 2}}), zeros());
		require(hi.Condition() > 1e23, "N8 high-conditioned fixture failed");
		emit("N8_C", hi.Values());
		emit("N8_E", hi.NumericalError());
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				require(hi.NumericalError()[i][j] < (i == 0 && j == 0 ? 1e12 : 1e-12), "N8 absolute response uncertainty uncontrolled");
		auto e = zeros();
		e[0][0] = 1e-24;
		refuse([&]
			   { ChargeNeutralNumberSusceptibility::Compute(full({{1e-24, 0, 0}, {0, 1, 0}, {0, 0, 2}}), e); }, "SupportedModeInsideUncertainty");
		refuse([&]
			   { ChargeNeutralNumberSusceptibility::Compute(full({{-1, 0, 0}, {0, 1, 0}, {0, 0, 2}}), zeros()); }, "Indefinite");
		refuse([&]
			   { ChargeNeutralNumberSusceptibility::Compute(full({{1, .1, 0}, {0, 1, 0}, {0, 0, 2}}), zeros()); }, "Nonsymmetric");
		TrackRFreeGasThermodynamicProvider p;
		for (double n : {0., p.NeutronOnsetBaryonDensityFm3(), p.MuonOnsetBaryonDensityFm3()})
			refuse([&]
				   { ChargeNeutralNumberSusceptibility::Compute(p.EquilibriumAt(n), {}); }, "ValueOnlyBoundary");
		for (auto nd : {std::pair<double, int>{.5, 3}, {.1, 2}, {1e-10, 1}})
		{
			auto a = ChargeNeutralNumberSusceptibility::Compute(p.EquilibriumAt(nd.first), zeros(nd.second));
			require(a.Support().size() == std::size_t(nd.second), "active dimension");
		}
		// A two-ended square-root susceptibility has exact integral pi/8 after
		// cancelling the spherical weight with the declared synthetic local C.
		auto r = request();
		r.intervals = {{0, 1, true, true, false}};
		r.local = [](double x)
		{M c=zeros();double f=std::sqrt(x*(1-x))/(4*M_PI*x*x)*1e-54;for(int i=0;i<3;++i)c[i][i]=f;return ChemicalLocalResponse{c,zeros(),{NumberAxis::Neutron,NumberAxis::Electron,NumberAxis::Muon}}; };
		auto both = GlobalChemicalNumberResponse::Compute(r);
		require(std::abs(both.Values()[0][0] - M_PI / 8) < 1e-12, "N2 measure/support");
		require(both.Partition().size() == 3 && both.Partition()[1] == .5, "N2 midpoint missing");
		r.intervals[0].first_order_discontinuity = true;
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "FirstOrderInterface");
		// Full refusal interval is excluded; a local evaluator that is queried there
		// throws, exposing either missing edge. Analytic omitted response = its volume.
		r = request();
		r.intervals = {{0, .25}, {.25, .75}, {.75, 1}};
		r.local = [](double x)
		{require(x<.4||x>.6,"N9 excluded node queried");M c=zeros();for(int i=0;i<3;++i)c[i][i]=1e-54;return ChemicalLocalResponse{c,zeros(),{NumberAxis::Neutron,NumberAxis::Electron,NumberAxis::Muon}}; };
		ChemicalRefusalCertificate f{.4, .6, .25, .75, 1, 1, "synthetic interior", "analytic unavailable interval", "flat exact geometry", .75, 0, 0, zeros()};
		double omitted = 4 * M_PI / 3 * (std::pow(.6, 3) - std::pow(.4, 3));
		for (int i = 0; i < 3; ++i)
			f.numerical_error[i][i] = omitted;
		r.refusals = {f};
		r.absolute_goal = M(3, std::vector<double>(3, 1));
		auto excluded = GlobalChemicalNumberResponse::Compute(r);
		require(std::abs(excluded.Values()[0][0] + omitted - 4 * M_PI / 3) < 1e-12, "N9 lost measure");
		r.refusals[0].containing_right = .5;
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "NotContained");
		r.refusals[0] = f;
		r.refusals[0].first_cell = 2;
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "NotContained");
		r.refusals[0] = f;
		r.refusals[0].availability_authority.clear();
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "AuthorityUnavailable");
		r.refusals[0] = f;
		r.refusals[0].nu_lower = 1;
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "GeometryExtremaContradicted");
		r = request();
		r.background = [](double x)
		{ return ChemicalMetric{.1 * x * x * x, 0}; };
		r.absolute_goal = M(3, std::vector<double>(3, 1));
		r.tail_error[1][1] = .1;
		r.pe_tail = ChemicalPeTailCertificate{1, .1, 0, 1e-4, 1e-8, 1.002, .1000001, 1.001, 1e-54, 1e-5, .1, "analytic positive pe comparison fixture"};
		(void)GlobalChemicalNumberResponse::Compute(r);
		auto saved_tail = *r.pe_tail;
		r.pe_tail->radius_upper = 1;
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "RadiusNotEnclosed");
		r.pe_tail = saved_tail;
		r.pe_tail->total_mass_upper = .1;
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "MassNotEnclosed");
		r.pe_tail = saved_tail;
		r.pe_tail->cut_mass = .09;
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "CutMismatch");
		r = request();
		auto g = std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(r));
		auto z = ChemicalImbalanceResponse::Compute(g, M(2, std::vector<double>(2, 1e-8)), M(2, std::vector<double>(2, 1e-8)));
		emit("GC11_Z", z.Values());
		// Independent rational two-zone fixture: input G_x values converted to y
		// by explicit integer algebra in this test; global Schur stays production.
		auto two = request(2);
		two.intervals = {{0, 1}, {1, 2}};
		two.local = [](double x)
		{M gx=x<1?M{{3,1,0},{1,2,0},{0,0,1}}:M{{2,0,1},{0,1,0},{1,0,3}};
  M t{{1,-1,-1},{0,1,0},{0,0,1}},cy=zeros();for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int a=0;a<3;++a)for(int b=0;b<3;++b)cy[i][j]+=t[i][a]*gx[a][b]*t[j][b];
  for(auto &row:cy)for(auto &v:row)v*=1e-54/(4*M_PI*x*x);
  return ChemicalLocalResponse{cy,zeros(),{NumberAxis::Neutron,NumberAxis::Electron,NumberAxis::Muon}}; };
		auto tg = std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(two));
		auto tz = ChemicalImbalanceResponse::Compute(tg, M(2, std::vector<double>(2, 1e-8)), M(2, std::vector<double>(2, 1e-8)));
		emit("GC10_Q", tz.ReductionDiagnostic().q);
		emit("GC10_Z", tz.Values());
		auto source = request();
		source.local = [](double x)
		{M c{{68,-5,-11},{-5,51,-13},{-11,-13,34}};for(auto &row:c)for(auto &v:row)v*=1e-54/(313*4*M_PI*x*x);return ChemicalLocalResponse{c,zeros(),{NumberAxis::Neutron,NumberAxis::Electron,NumberAxis::Muon}}; };
		auto sg = std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(source));
		auto sz = ChemicalImbalanceResponse::Compute(sg, M(2, std::vector<double>(2, 1e-8)), M(2, std::vector<double>(2, 1e-8)));
		emit("GC11_SOURCE_Z", sz.Values());
		// Physical support yields one beta channel for npe, and none for all-pe.
		auto lower = request();
		lower.local = [](double)
		{M c=zeros();c[0][0]=c[1][1]=1e-54;return ChemicalLocalResponse{c,zeros(),{NumberAxis::Neutron,NumberAxis::Electron}}; };
		auto lg = std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(lower));
		auto lz = ChemicalImbalanceResponse::Compute(lg, {{1e-8}}, {{1e-8}});
		require(lz.Channels() == std::vector<ImbalanceChannel>{ImbalanceChannel::Npe}, "npe channel set");
		lower.local = [](double)
		{M c=zeros();c[1][1]=1e-54;return ChemicalLocalResponse{c,zeros(),{NumberAxis::Electron}}; };
		auto pg = std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(lower));
		refuse([&]
			   { ChemicalImbalanceResponse::Compute(pg, {{1e-8}}, {{1e-8}}); }, "NoSupportedBetaChannel");
		auto saved = *r.lifetime->revision;
		for (auto member : {&ChemicalRevision::model, &ChemicalRevision::provider_revision, &ChemicalRevision::provider_bytes, &ChemicalRevision::background, &ChemicalRevision::domain, &ChemicalRevision::partition, &ChemicalRevision::onset, &ChemicalRevision::tail, &ChemicalRevision::accuracy, &ChemicalRevision::lifetime_token})
		{
			r.lifetime->revision.get()->*member += "changed";
			refuse([&]
				   { z.Values(); }, "StaleChemical");
			*r.lifetime->revision = saved;
		}
		r.lifetime->revision->alive = false;
		refuse([&]
			   { z.Values(); }, "ExpiredLifetime");
		r.lifetime->revision->alive = true;
		r.lifetime.reset();
		g.reset();
		(void)z.Values();
		r = request();
		for (int i = 0; i < 3; ++i)
			r.background_error[i][i] = .11 * (4 * M_PI / 3);
		r.absolute_goal = M(3, std::vector<double>(3, 1));
		auto unresolved_q = std::make_shared<GlobalChemicalNumberResponse>(GlobalChemicalNumberResponse::Compute(r));
		refuse([&]
			   { ChemicalImbalanceResponse::Compute(unresolved_q, M(2, std::vector<double>(2, 10)), M(2, std::vector<double>(2, 10))); }, "SupportedModeInsideUncertainty");
		r = request();
		r.background = [](double x)
		{ return ChemicalMetric{.5 * x * (1 - std::numeric_limits<double>::epsilon()), 0}; };
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "GeometryArithmeticUnresolved");
		r = request();
		r.local = [](double)
		{ return ChemicalLocalResponse{zeros(), zeros(), {static_cast<NumberAxis>(99)}}; };
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "InvalidPhysicalSupportAxis");
		r = request();
		r.structural_zeros = {{NumberAxis::Neutron, NumberAxis::Electron}};
		r.structural_zero_authority = "analytic diagonal susceptibility";
		require(GlobalChemicalNumberResponse::Compute(r).Diagnostics().structural_zeros == r.structural_zeros, "structural-zero declaration lost");
		r.structural_zeros = {{NumberAxis::Neutron, NumberAxis::Neutron}};
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "StructuralZeroContradicted");
		r = request();
		r.background_error = M(3, std::vector<double>(3, 10));
		r.absolute_goal = r.background_error;
		refuse([&]
			   { GlobalChemicalNumberResponse::Compute(r); }, "SupportedModeInsideUncertainty");
		std::cout << "CONTRACTS PASS\n";
		return 0;
	}
	catch (const std::exception &e)
	{
		std::cerr << "FAIL " << e.what() << '\n';
		return 1;
	}
}
