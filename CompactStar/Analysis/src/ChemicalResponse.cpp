#include <CompactStar/Analysis/ChemicalResponse.hpp>
#include <CompactStar/AngularVelocity.hpp>
#include <CompactStar/Geometry.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace CompactStar::Analysis
{
namespace
{
using M = ChemicalMatrix;
constexpr double u = std::numeric_limits<double>::epsilon();
void need(bool ok, const char *why)
{
	if (!ok)
		throw std::runtime_error(why);
}
M zero(std::size_t n, std::size_t m) { return M(n, std::vector<double>(m)); }
M eye(std::size_t n)
{
	auto a = zero(n, n);
	for (std::size_t i = 0; i < n; ++i)
		a[i][i] = 1;
	return a;
}
void shape(const M &a, std::size_t n, std::size_t m, bool positive = false)
{
	need(a.size() == n, "ChemicalMatrixRowMismatch");
	for (const auto &row : a)
	{
		need(row.size() == m, "ChemicalMatrixColumnMismatch");
		for (double x : row)
			need(std::isfinite(x) && (!positive || x >= 0), "InvalidChemicalMatrixEntry");
	}
}
M absM(M a)
{
	for (auto &row : a)
		for (auto &v : row)
			v = std::abs(v);
	return a;
}
M trans(const M &a)
{
	auto b = zero(a[0].size(), a.size());
	for (std::size_t i = 0; i < a.size(); ++i)
		for (std::size_t j = 0; j < a[0].size(); ++j)
			b[j][i] = a[i][j];
	return b;
}
M plus(M a, const M &b, double scale = 1)
{
	shape(b, a.size(), a[0].size());
	for (std::size_t i = 0; i < a.size(); ++i)
		for (std::size_t j = 0; j < a[0].size(); ++j)
			a[i][j] += scale * b[i][j];
	return a;
}
M times(M a, double s)
{
	for (auto &row : a)
		for (auto &v : row)
			v *= s;
	return a;
}
struct Sum
{
	double s = 0, c = 0;
	void add(double v)
	{
		double t = s + v;
		c += std::abs(s) >= std::abs(v) ? (s - t) + v : (v - t) + s;
		s = t;
	}
	double value() const { return s + c; }
};
M mul(const M &a, const M &b)
{
	need(!a.empty() && !b.empty() && a[0].size() == b.size(), "ChemicalProductShape");
	auto c = zero(a.size(), b[0].size());
	for (std::size_t i = 0; i < c.size(); ++i)
		for (std::size_t j = 0; j < c[0].size(); ++j)
		{
			Sum s;
			for (std::size_t k = 0; k < b.size(); ++k)
				s.add(a[i][k] * b[k][j]);
			c[i][j] = s.value();
		}
	return c;
}
double norm(const M &a)
{
	double v = 0;
	for (const auto &r : a)
	{
		double s = 0;
		for (double x : r)
			s += std::abs(x);
		v = std::max(v, s);
	}
	return v;
}
std::vector<double> spectrum(const M &a)
{
	const auto n = a.size();
	auto g = gsl_matrix_alloc(n, n);
	auto e = gsl_vector_alloc(n);
	auto w = gsl_eigen_symm_alloc(n);
	need(g && e && w, "EigenAllocationFailure");
	for (std::size_t i = 0; i < n; ++i)
		for (std::size_t j = 0; j < n; ++j)
			gsl_matrix_set(g, i, j, (a[i][j] + a[j][i]) / 2);
	const int rc = gsl_eigen_symm(g, e, w);
	std::vector<double> v(n);
	for (std::size_t i = 0; i < n; ++i)
		v[i] = gsl_vector_get(e, i);
	gsl_eigen_symm_free(w);
	gsl_vector_free(e);
	gsl_matrix_free(g);
	need(rc == 0, "EigenSolveFailure");
	std::sort(v.begin(), v.end());
	return v;
}
// Pivoted elimination is used only for the uncertainty majorant I-R.
// Physical symmetric positive response matrices use scaled Cholesky below.
M general_solve(M a, M b)
{
	const auto n = a.size();
	for (std::size_t k = 0; k < n; ++k)
	{
		std::size_t p = k;
		for (std::size_t j = k + 1; j < n; ++j)
			if (std::abs(a[j][k]) > std::abs(a[p][k]))
				p = j;
		need(a[p][k] != 0, "UnresolvedUncertaintyMajorant");
		std::swap(a[p], a[k]);
		std::swap(b[p], b[k]);
		const double d = a[k][k];
		for (std::size_t j = k; j < n; ++j)
			a[k][j] /= d;
		for (auto &v : b[k])
			v /= d;
		for (std::size_t i = 0; i < n; ++i)
			if (i != k)
			{
				const double f = a[i][k];
				for (std::size_t j = k; j < n; ++j)
					a[i][j] -= f * a[k][j];
				for (std::size_t j = 0; j < b[i].size(); ++j)
					b[i][j] -= f * b[k][j];
			}
	}
	return b;
}
struct Inverse
{
	M x, error, solve_error;
	double rho, condition;
};
Inverse inverse_spd(M a, const M &e)
{
	const auto n = a.size();
	need(n > 0 && n <= 3, "UnsupportedActiveDimension");
	shape(a, n, n);
	shape(e, n, n, true);
	std::vector<double> d(n);
	M scaled = a, se = e;
	for (std::size_t i = 0; i < n; ++i)
	{
		need(a[i][i] > 0, "IndefiniteChemicalResponse");
		d[i] = std::sqrt(a[i][i]);
	}
	for (std::size_t i = 0; i < n; ++i)
		for (std::size_t j = 0; j < n; ++j)
		{
			need(std::abs(a[i][j] - a[j][i]) <= e[i][j] + e[j][i] + 32 * u * (std::abs(a[i][j]) + std::abs(a[j][i])), "NonsymmetricBeyondBudget");
			scaled[i][j] = (a[i][j] + a[j][i]) / 2 / d[i] / d[j];
			se[i][j] = e[i][j] / d[i] / d[j] + 32 * u * std::abs(scaled[i][j]);
		}
	const auto ev = spectrum(scaled);
	need(ev.front() > norm(se), "SupportedModeInsideUncertainty");
	auto l = zero(n, n);
	for (std::size_t i = 0; i < n; ++i)
		for (std::size_t j = 0; j <= i; ++j)
		{
			Sum s;
			s.add(scaled[i][j]);
			for (std::size_t k = 0; k < j; ++k)
				s.add(-l[i][k] * l[j][k]);
			if (i == j)
			{
				need(s.value() > 0, "IndefiniteChemicalResponse");
				l[i][j] = std::sqrt(s.value());
			}
			else
				l[i][j] = s.value() / l[j][j];
		}
	auto x = zero(n, n);
	for (std::size_t col = 0; col < n; ++col)
	{
		std::vector<double> y(n), z(n);
		for (std::size_t i = 0; i < n; ++i)
		{
			Sum s;
			s.add(i == col ? 1 : 0);
			for (std::size_t j = 0; j < i; ++j)
				s.add(-l[i][j] * y[j]);
			y[i] = s.value() / l[i][i];
		}
		for (int i = int(n) - 1; i >= 0; --i)
		{
			Sum s;
			s.add(y[i]);
			for (std::size_t j = i + 1; j < n; ++j)
				s.add(-l[j][i] * z[j]);
			z[i] = s.value() / l[i][i];
		}
		for (std::size_t i = 0; i < n; ++i)
			x[i][col] = z[i] / d[i] / d[col];
	}
	// The computed solve residual, dot-product arithmetic and uncertainty in H
	// all enter the componentwise inverse perturbation majorant explicitly.
	auto ax = absM(x), aa = absM(a);
	// Left residual: X(A+dA)=I-D+X*dA. Retain the complete
	// Neumann majorant, including solve residual/input-error cross terms.
	auto residual = plus(absM(plus(eye(n), mul(x, a), -1)), times(mul(ax, aa), 32 * u));
	need(norm(residual) < 1, "GlobalSolveResidualUnresolved");
	auto solve_error = general_solve(plus(eye(n), residual, -1), mul(residual, ax));
	auto r = plus(mul(ax, e), residual);
	const double rho = norm(r);
	need(rho < 1, "InversePerturbationUnresolved");
	auto error = general_solve(plus(eye(n), r, -1), mul(r, ax));
	for (auto &row : error)
		for (auto &v : row)
		{
			need(v >= 0 && std::isfinite(v), "InvalidInverseErrorMajorant");
			v = std::nextafter(v * (1 + 64 * u), std::numeric_limits<double>::infinity());
		}
	return {x, error, solve_error, rho, norm(a) * norm(x)};
}
void goal(const M &e, const M &g)
{
	shape(g, e.size(), e[0].size(), true);
	for (std::size_t i = 0; i < e.size(); ++i)
		for (std::size_t j = 0; j < e[i].size(); ++j)
			need(e[i][j] <= g[i][j], "AccuracyGoalUnmet");
}
std::string bytes(const std::string &p)
{
	std::ifstream in(p, std::ios::binary);
	need(bool(in), "SourceBytesUnavailable");
	std::ostringstream contents;
	contents << in.rdbuf();
	need(!in.bad(), "SourceReadFailure");
	return contents.str();
}
std::string provider_text(const ILocalThermodynamicProvider &p)
{
	const auto &m = p.Metadata();
	return m.model_id + '\n' + m.model_revision + '\n' + m.particle_content + '\n' + m.coordinate_chart + '\n' + m.temperature_scope + '\n' + m.rest_mass_convention + '\n' + m.lepton_ownership + '\n' + m.smooth_domain;
}
std::string revision(const ChemicalLifetime &l)
{
	need(l.revision && l.revision->alive, "ExpiredLifetimeToken");
	return l.revision->Serialize() + (l.provider ? provider_text(*l.provider) : "injected local evaluator");
}
M supported(const M &a, const std::vector<NumberAxis> &s)
{
	auto b = zero(s.size(), s.size());
	for (std::size_t i = 0; i < s.size(); ++i)
		for (std::size_t j = 0; j < s.size(); ++j)
			b[i][j] = a[int(s[i])][int(s[j])];
	return b;
}
ChemicalReductionDiagnostic reduce(const GlobalChemicalNumberResponse &g)
{
	const auto &s = g.Support();
	need(!s.empty() && s[0] == NumberAxis::Neutron, "NoSupportedBetaChannel");
	need(s.size() >= 2, "NoSupportedLeptonChannel");
	auto a = supported(g.Values(), s), e = supported(g.NumericalError(), s);
	auto U = eye(s.size());
	for (auto &v : U[0])
		v = 1;
	auto gx = mul(U, mul(a, trans(U))), ex = mul(U, mul(e, trans(U)));
	ex = plus(ex, times(mul(absM(U), mul(absM(a), trans(absM(U)))), 32 * u));
	const double den = gx[0][0], de = ex[0][0];
	need(den > de, "GlobalBaryonDenominatorUnresolved");
	const std::size_t n = s.size() - 1;
	auto q = zero(n, n), qe = q, ar = q;
	for (std::size_t i = 0; i < n; ++i)
		for (std::size_t j = 0; j < n; ++j)
		{
			const double h = gx[i + 1][0], k = gx[0][j + 1], eh = ex[i + 1][0], ek = ex[0][j + 1];
			q[i][j] = gx[i + 1][j + 1] - h * k / den;
			ar[i][j] = 32 * u * (std::abs(gx[i + 1][j + 1]) + std::abs(h * k / den));
			qe[i][j] = ex[i + 1][j + 1] + (std::abs(h) * ek + eh * std::abs(k) + eh * ek) / (den - de) + std::abs(h * k) * de / (den * (den - de)) + ar[i][j];
		}
	(void)inverse_spd(q, qe);
	return {q, qe, ar, den, de};
}
std::string structural_text(const FixedBaryonNumberResponse &s)
{
	std::ostringstream o;
	// Length prefixes prevent identity/number concatenation collisions.
	auto item = [&](const auto &value)
	{std::ostringstream part;part<<std::setprecision(17)<<value;auto text=part.str();o<<text.size()<<':'<<text; };
	const auto &m = s.metadata;
	item(m.domain.boundary_definition);
	item(int(m.domain.type));
	item(m.domain.inner_pressure_km_minus2);
	item(m.domain.outer_pressure_km_minus2);
	item(m.central_state_definition);
	item(m.central_energy_km_minus2);
	item(m.sequence_radial_resolution);
	item(m.sequence_branch_policy);
	item(m.tail_callback_identity);
	item(m.tail_callback_revision);
	item(m.refinement_policy);
	item(m.q_normalization);
	item(m.profile_units);
	item(m.differentiation_coordinate);
	item(m.surface_authority);
	item(int(m.surface));
	for (double value : {s.A_B, s.B_B, s.B_B_error, s.central_energy_per_q, s.central_energy_per_q_error, s.conditioning, s.baryon_residual, s.baryon_budget, s.charge_residual, s.charge_budget})
		item(value);
	item(m.sources.size());
	for (const auto &p : m.sources)
	{
		item(p.eos_snapshot.identity);
		item(p.eos_snapshot.revision);
		item(p.eos_snapshot.physical_domain);
		item(p.table_contents);
		item(p.profile_version);
	}
	item(m.species.size());
	for (const auto &p : m.species)
	{
		item(p.label);
		item(p.baryon_number);
		item(p.charge);
	}
	item(m.tail_inputs.size());
	for (const auto &t : m.tail_inputs)
	{
		item(int(t.semantics));
		item(t.authority);
		item(t.policy_revision);
		for (const auto *v : {&t.count_correction, &t.count_error, &t.response_correction, &t.response_error})
		{
			item(v->size());
			for (double x : *v)
				item(x);
		}
	}
	item(m.profile_node_counts.size());
	for (auto x : m.profile_node_counts)
		item(x);
	for (const auto *v : {&m.steps, &m.achieved_central_energy_km_minus2})
	{
		item(v->size());
		for (double x : *v)
			item(x);
	}
	const auto &values = s.Values(), &errors = s.Errors();
	item(values.size());
	for (double x : values)
		item(x);
	item(errors.size());
	for (double x : errors)
		item(x);
	return o.str();
}
void coverage(const ChemicalLifetime &l, const FixedBaryonNumberResponse &s)
{
	need(s.metadata.sources.size() >= 16, "IncompleteStructuralDependencyBundle");
	for (const auto &p : s.metadata.sources)
	{
		auto it = std::find_if(l.stars.begin(), l.stars.end(), [&](const auto &v)
							   { return v && v.get() == p.star; });
		need(it != l.stars.end(), "ForeignOrUnownedStructuralSource");
	}
}
} // namespace

ChargeNeutralNumberSusceptibility ChargeNeutralNumberSusceptibility::Compute(const ActiveLocalThermodynamicEvaluation &a, const M &he)
{
	ChargeNeutralNumberSusceptibility out;
	std::visit([&](const auto &v)
			   {using V=std::decay_t<decltype(v)>;constexpr auto n=V::response_dimension;
        if constexpr(n==0) throw std::runtime_error("ValueOnlyBoundaryHasNoHessian");
        else {M h=zero(n,n);for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)h[i][j]=v.hessian(i,j);
            // Solve in the canonical active y chart after an integer congruence.
            // This is algebraically T H_x^-1 T^t, with no inverse cancellation
            // of the exactly decoupled neutron/lepton free-gas blocks. Each
            // performed addition contributes its own absolute rounding term.
            shape(he,n,n,true);auto U=eye(n);if constexpr(n>1)for(auto &x:U[0])x=1;
            M hy=zero(n,n),ey=zero(n,n);
            for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){double value=0,ar=0;
                for(std::size_t a=0;a<n;++a)for(std::size_t b=0;b<n;++b)if(U[a][i]&&U[b][j]){
                    const double term=h[a][b];const double next=value+term;
                    // Addition of zero and exact opposite binary values is exact.
                    if(value!=0&&term!=0&&value!=-term)ar+=u*std::abs(next)/(1-u);
                    value=next;ey[i][j]+=he[a][b];}
                hy[i][j]=value;ey[i][j]+=ar;}
            auto inv=inverse_spd(hy,ey);M t=zero(3,n);
            if constexpr(n==3){t={{1,0,0},{0,1,0},{0,0,1}};out.support_={NumberAxis::Neutron,NumberAxis::Electron,NumberAxis::Muon};}
            if constexpr(n==2){t={{1,0},{0,1},{0,0}};out.support_={NumberAxis::Neutron,NumberAxis::Electron};}
            if constexpr(n==1){t={{0},{1},{0}};out.support_={NumberAxis::Electron};}
            out.c_=mul(t,mul(inv.x,trans(t)));out.error_=mul(t,mul(inv.error,trans(t)));
            out.rho_=inv.rho;out.condition_=inv.condition;
        } }, a);
	return out;
}
std::string ChemicalRevision::Serialize() const
{
	std::ostringstream o;
	for (const auto *s : {&model, &provider_revision, &provider_bytes, &particle_constants, &background, &metric, &basis, &domain, &partition, &onset, &tail, &accuracy, &lifetime_token})
	{
		need(!s->empty(), "MissingChemicalProvenance");
		o << s->size() << ':' << *s;
	}
	return o.str();
}

GlobalChemicalNumberResponse GlobalChemicalNumberResponse::Compute(const ChemicalIntegrationRequest &r)
{
	need(r.lifetime && r.background && r.local, "MissingIntegrationDependency");
	need(!r.intervals.empty(), "EmptyChemicalPartition");
	need(r.order == 8 || r.order == 16 || r.order == 32 || r.order == 64, "UnsupportedGaussOrder");
	for (const M *e : {&r.center_error, &r.tail_error, &r.background_error, &r.absolute_goal})
		shape(*e, 3, 3, true);
	need(!r.center_authority.empty() && !r.tail_authority.empty() && !r.background_error_authority.empty(), "MissingCenterTailBackgroundAuthority");
	if (r.pe_tail)
	{
		const auto &t = *r.pe_tail;
		for (double v : {t.cut_radius, t.cut_mass, t.h_upper, t.epsilon_upper, t.bootstrap_radius_upper, t.total_mass_upper, t.radius_upper, t.susceptibility_upper, t.lapse_error_upper, t.response_upper})
			need(std::isfinite(v) && v > 0, "InvalidPeTailCertificate");
		need(std::isfinite(t.cut_nu) && !t.source_hypotheses.empty(), "MissingPeTailSourceAuthority");
		const long double R = t.cut_radius, M = t.cut_mass, H = t.h_upper, mu = t.total_mass_upper;
		const auto cut_metric = r.background(t.cut_radius);
		need(t.cut_radius == r.intervals.back().right && cut_metric.mass_km == t.cut_mass && cut_metric.nu == t.cut_nu, "PeTailCutMismatch");
		const long double pi = std::acos(-1.L), den = 2 * M / R - (1 - 2 * M / R) * std::expm1(2 * H);
		need(R > 2 * mu && 1 - H * R / M > 0 && den > 0, "PeTailGeometryUnresolved");
		need(t.bootstrap_radius_upper >= R / (1 - H * R / M), "PeTailBootstrapNotEnclosed");
		const long double rb = t.bootstrap_radius_upper;
		need(mu >= M + 4 * pi / 3 * t.epsilon_upper * ((rb - R) * (rb * rb + rb * R + R * R)), "PeTailMassNotEnclosed");
		// Known cut mass remains in the boundary factor; total mass UPPER is
		// mandatory in the numerator. This check rejects the old-M mutation.
		need(t.radius_upper >= 2 * mu / den && t.radius_upper <= rb, "PeTailRadiusNotEnclosed");
		const long double ru = t.radius_upper;
		need(t.lapse_error_upper >= (mu - M) / (R * (1 - 2 * mu / R)), "PeTailLapseNotEnclosed");
		const long double bound = 4 * pi / 3 * ((ru - R) * (ru * ru + ru * R + R * R)) * std::exp(-static_cast<long double>(t.cut_nu) + t.lapse_error_upper) * t.susceptibility_upper / std::sqrt(1 - 2 * mu / R) * 1e54L;
		need(t.response_upper >= bound && r.tail_error[1][1] >= t.response_upper, "PeTailResponseNotEnclosed");
	}
	need(r.structural_zeros.empty() || !r.structural_zero_authority.empty(), "StructuralZeroAuthorityUnavailable");
	for (const auto &entry : r.structural_zeros)
		need(int(entry.first) >= 0 && int(entry.first) < 3 && int(entry.second) >= 0 && int(entry.second) < 3, "InvalidStructuralZeroAxis");
	GlobalChemicalNumberResponse out;
	out.lifetime_ = r.lifetime;
	out.snapshot_ = revision(*r.lifetime);
	for (const auto &p : r.lifetime->source_files)
		out.file_contents_.push_back(bytes(p));
	for (const auto &s : r.lifetime->stars)
	{
		need(bool(s), "NullStarOwner");
		NumberProvenance p;
		p.star = s.get();
		p.profile = &s->Profile();
		p.profile_version = s->Profile().Version();
		out.profiles_.push_back(p);
	}
	std::vector<ChemicalInterval> partition;
	for (std::size_t i = 0; i < r.intervals.size(); ++i)
	{
		const auto &s = r.intervals[i];
		need(s.left >= 0 && s.right > s.left, "InvalidChemicalInterval");
		need(!s.first_order_discontinuity, "FirstOrderInterfaceLawUnavailable");
		if (i)
			need(r.intervals[i - 1].right == s.left, "UnrepresentedChemicalPartitionGap");
		std::vector<double> cuts{s.left, s.right};
		if (s.left_continuous_onset && s.right_continuous_onset)
			cuts.push_back(s.left + (s.right - s.left) / 2);
		for (const auto &f : r.refusals)
			for (double edge : {f.left, f.right})
				if (edge > s.left && edge < s.right)
					cuts.push_back(edge);
		std::sort(cuts.begin(), cuts.end());
		cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
		for (std::size_t j = 1; j < cuts.size(); ++j)
			partition.push_back({cuts[j - 1], cuts[j], cuts[j - 1] == s.left && s.left_continuous_onset, cuts[j] == s.right && s.right_continuous_onset, false});
	}
	M refusal_error = zero(3, 3);
	for (const auto &f : r.refusals)
	{
		need(f.right > f.left && f.first_cell <= f.last_cell && f.containing_left <= f.left && f.right <= f.containing_right, "RefusalIntervalNotContained");
		need(f.last_cell < r.intervals.size() && f.containing_left == r.intervals[f.first_cell].left && f.containing_right == r.intervals[f.last_cell].right, "RefusalContainingCellsUnauthenticated");
		need(f.left >= r.intervals.front().left && f.right <= r.intervals.back().right, "RefusalOutsideDomain");
		need(!f.branch.empty() && !f.availability_authority.empty() && !f.geometry_extrema_authority.empty(), "RefusalAuthorityUnavailable");
		need(f.radius_upper >= f.right && f.mass_upper >= 0 && f.left > 2 * f.mass_upper && std::isfinite(f.nu_lower), "RefusalGeometryUnresolved");
		for (std::size_t cell = f.first_cell; cell <= f.last_cell; ++cell)
			for (double edge : {r.intervals[cell].left, r.intervals[cell].right})
			{
				auto metric = r.background(edge);
				need(edge <= f.radius_upper && metric.mass_km <= f.mass_upper && metric.nu >= f.nu_lower, "RefusalGeometryExtremaContradicted");
			}
		shape(f.numerical_error, 3, 3, true);
		refusal_error = plus(refusal_error, f.numerical_error);
	}
	for (const auto &s : partition)
		out.partition_.push_back(s.left);
	out.partition_.push_back(partition.back().right);
	struct Integral
	{
		M value, error;
		std::array<bool, 3> support{};
		std::size_t nodes = 0;
	};
	auto integrate = [&](unsigned order, bool bisect)
	{
		Integral ans{zero(3, 3), zero(3, 3)};
		std::array<Sum, 9> sums, errors;
		auto *rule = gsl_integration_glfixed_table_alloc(order);
		need(rule, "GaussRuleUnavailable");
		try
		{
			for (const auto &s : partition)
			{
				bool excluded = false;
				for (const auto &f : r.refusals)
					if (s.left >= f.left && s.right <= f.right)
						excluded = true;
				if (excluded)
					continue;
				for (unsigned part = 0; part < (bisect ? 2u : 1u); ++part)
				{
					double a = s.left, b = s.right;
					if (bisect)
					{
						double mid = a + (b - a) / 2;
						if (part)
							b = s.right, a = mid;
						else
							b = mid;
					}
					const bool left = a == s.left && s.left_continuous_onset, right = b == s.right && s.right_continuous_onset;
					need(!(left && right), "DoubleOnsetNotSplit");
					for (unsigned j = 0; j < order; ++j)
					{
						double t, w;
						gsl_integration_glfixed_point(0, 1, j, &t, &w, rule);
						double x = a + (b - a) * t, jac = b - a;
						if (left)
						{
							x = a + (b - a) * t * t;
							jac *= 2 * t;
						}
						if (right)
						{
							x = b - (b - a) * t * t;
							jac *= 2 * t;
						}
						need(x > a && x < b, "MappedNodeNotRepresentableInsideSegment");
						auto metric = r.background(x);
						const double compactness = 2 * metric.mass_km / x, geometry_denominator = 1 - compactness;
						const double geometry_rounding = 4 * u * std::abs(compactness);
						need(std::isfinite(metric.nu) && std::isfinite(compactness) && geometry_denominator > geometry_rounding, "GeometryArithmeticUnresolved");
						// Arithmetic in 1-2m/r is amplified near a horizon. It
						// must not be hidden inside a uniform relative floor.
						const double measure_arithmetic = 64 * u + geometry_rounding / (geometry_denominator - geometry_rounding);
						const double measure = 1e54 * Geometry::ProperVolumeWeight(x, metric.mass_km) * std::exp(-metric.nu) * jac * w;
						need(std::isfinite(measure) && measure > 0, "InvalidChemicalMeasure");
						auto local = r.local(x);
						shape(local.c, 3, 3);
						shape(local.numerical_error, 3, 3, true);
						for (const auto &entry : r.structural_zeros)
							need(local.c[int(entry.first)][int(entry.second)] == 0, "StructuralZeroContradicted");
						std::array<bool, 3> support{};
						for (auto axis : local.support)
						{
							need(int(axis) >= 0 && int(axis) < 3, "InvalidPhysicalSupportAxis");
							support[int(axis)] = true;
							ans.support[int(axis)] = true;
						}
						for (int i = 0; i < 3; ++i)
							for (int k = 0; k < 3; ++k)
							{
								if (!support[i] || !support[k])
									need(local.c[i][k] == 0 && local.numerical_error[i][k] == 0, "FabricatedInactiveSupport");
								sums[3 * i + k].add(measure * local.c[i][k]);
								errors[3 * i + k].add(measure * (local.numerical_error[i][k] + measure_arithmetic * std::abs(local.c[i][k])));
							}
						++ans.nodes;
					}
				}
			}
		}
		catch (...)
		{
			gsl_integration_glfixed_table_free(rule);
			throw;
		}
		gsl_integration_glfixed_table_free(rule);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
			{
				ans.value[i][j] = sums[3 * i + j].value();
				ans.error[i][j] = errors[3 * i + j].value();
			}
		return ans;
	};
	const auto a = integrate(r.order / 2 < 8 ? 8 : r.order / 2, false), b = integrate(r.order, false), c = integrate(r.order * 2 <= 64 ? r.order * 2 : 64, false), d = integrate(r.order * 2 <= 64 ? r.order * 2 : 64, true);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
		{
			const double coarse = std::abs(a.value[i][j] - b.value[i][j]), fine = std::abs(b.value[i][j] - c.value[i][j]);
			const double floor = a.error[i][j] + b.error[i][j] + c.error[i][j];
			need(fine <= coarse || fine <= floor, "QuadratureNotContractingAboveNumericalFloor");
		}
	out.diagnostics_.order = r.order;
	out.diagnostics_.structural_zeros = r.structural_zeros;
	out.diagnostics_.structural_zero_authority = r.structural_zero_authority;
	out.diagnostics_.validation_orders = {r.order / 2 < 8 ? 8 : r.order / 2, r.order, r.order * 2 <= 64 ? r.order * 2 : 64, r.order * 2 <= 64 ? r.order * 2 : 64};
	out.diagnostics_.validation_integrals = {a.value, b.value, c.value, d.value};
	out.diagnostics_.local_arithmetic_error = b.error;
	out.diagnostics_.center_error = r.center_error;
	out.diagnostics_.tail_error = r.tail_error;
	out.diagnostics_.background_error = r.background_error;
	out.diagnostics_.refusal_error = refusal_error;
	out.diagnostics_.absolute_goal = r.absolute_goal;
	out.diagnostics_.source_intervals = r.intervals;
	out.diagnostics_.refusal_certificates = r.refusals;
	out.diagnostics_.pe_tail = r.pe_tail;
	out.diagnostics_.center_authority = r.center_authority;
	out.diagnostics_.tail_authority = r.tail_authority;
	out.diagnostics_.background_error_authority = r.background_error_authority;
	out.diagnostics_.validation_node_count = a.nodes + b.nodes + c.nodes + d.nodes;
	out.g_ = b.value;
	out.nodes_ = b.nodes;
	out.quadrature_error_ = times(plus(plus(absM(plus(a.value, b.value, -1)), absM(plus(b.value, c.value, -1))), absM(plus(c.value, d.value, -1))), 2);
	out.error_ = plus(plus(plus(plus(plus(b.error, out.quadrature_error_), r.center_error), r.tail_error), r.background_error), refusal_error);
	for (int i = 0; i < 3; ++i)
		if (b.support[i])
			out.support_.push_back(static_cast<NumberAxis>(i));
	need(!out.support_.empty(), "NoGlobalPhysicalSupport");
	auto supported_g = supported(out.g_, out.support_), supported_e = supported(out.error_, out.support_);
	auto inv = inverse_spd(supported_g, supported_e);
	out.condition_ = inv.condition;
	out.eigenvalues_ = spectrum(supported_g);
	goal(out.error_, r.absolute_goal);
	out.RequireCurrent();
	return out;
}
void GlobalChemicalNumberResponse::RequireCurrent() const
{
	need(lifetime_ && revision(*lifetime_) == snapshot_, "StaleChemicalDependency");
	need(lifetime_->stars.size() == profiles_.size() && lifetime_->source_files.size() == file_contents_.size(), "LifetimeOwnerCoverageChanged");
	for (std::size_t i = 0; i < profiles_.size(); ++i)
	{
		const auto &p = profiles_[i];
		const auto &s = lifetime_->stars[i];
		need(s && s.get() == p.star && &s->Profile() == p.profile && s->Profile().Version() == p.profile_version, "StaleChemicalProfile");
	}
	for (std::size_t i = 0; i < file_contents_.size(); ++i)
		need(bytes(lifetime_->source_files[i]) == file_contents_[i], "StaleChemicalSourceBytes");
}
const M &GlobalChemicalNumberResponse::Values() const
{
	RequireCurrent();
	return g_;
}
const M &GlobalChemicalNumberResponse::NumericalError() const
{
	RequireCurrent();
	return error_;
}
const std::vector<NumberAxis> &GlobalChemicalNumberResponse::Support() const
{
	RequireCurrent();
	return support_;
}
const std::vector<double> &GlobalChemicalNumberResponse::Eigenvalues() const
{
	RequireCurrent();
	return eigenvalues_;
}
ChemicalImbalanceResponse ChemicalImbalanceResponse::Compute(std::shared_ptr<const GlobalChemicalNumberResponse> g, const M &q_goal, const M &z_goal)
{
	need(bool(g), "MissingGlobalChemicalDependency");
	g->RequireCurrent();
	const auto q = reduce(*g);
	goal(q.numerical_error, q_goal);
	const auto &s = g->Support();
	const std::size_t n = s.size(), m = n - 1;
	auto a = supported(g->Values(), s), e = supported(g->NumericalError(), s);
	auto inv = inverse_spd(a, e);
	auto l = zero(n, m);
	for (std::size_t j = 0; j < m; ++j)
	{
		l[0][j] = -1;
		l[j + 1][j] = 1;
	}
	ChemicalImbalanceResponse out;
	out.global_ = std::move(g);
	out.q_goal_ = q_goal;
	out.z_goal_ = z_goal;
	out.z_ = mul(trans(l), mul(inv.x, l));
	out.arithmetic_error_ = times(mul(trans(absM(l)), mul(absM(inv.x), absM(l))), 32 * u);
	out.global_solve_error_ = mul(trans(absM(l)), mul(inv.solve_error, absM(l)));
	out.error_ = plus(mul(trans(absM(l)), mul(inv.error, absM(l))), out.arithmetic_error_);
	out.rho_ = inv.rho;
	out.condition_ = inv.condition;
	for (std::size_t i = 1; i < n; ++i)
		out.channels_.push_back(s[i] == NumberAxis::Electron ? ImbalanceChannel::Npe : ImbalanceChannel::NpMu);
	goal(out.error_, z_goal);
	out.RequireCurrent();
	return out;
}
void ChemicalImbalanceResponse::RequireCurrent() const
{
	need(bool(global_), "MissingGlobalChemicalDependency");
	global_->RequireCurrent();
}
ChemicalReductionDiagnostic ChemicalImbalanceResponse::ReductionDiagnostic() const
{
	RequireCurrent();
	return reduce(*global_);
}
double ChemicalImbalanceResponse::PaperZ(ImbalanceChannel row, ImbalanceChannel col) const
{
	RequireCurrent();
	auto i = std::find(channels_.begin(), channels_.end(), row), j = std::find(channels_.begin(), channels_.end(), col);
	need(i != channels_.end() && j != channels_.end(), "UnsupportedPaperChannel");
	return z_[i - channels_.begin()][j - channels_.begin()];
}
RotochemicalSpinDrive RotochemicalSpinDrive::Compute(std::shared_ptr<const ChemicalImbalanceResponse> z, std::shared_ptr<const FixedBaryonNumberResponse> k, std::vector<double> vi, std::vector<double> gn, std::vector<double> gv)
{
	need(z && k, "MissingSpinDriveDependency");
	z->RequireCurrent();
	auto owner = z->Global()->Lifetime();
	coverage(*owner, *k);
	k->RequireCurrent();
	need(!owner->stars.empty() && k->metadata.sources.front().star == owner->stars.front().get(), "ForeignCentralStructuralSource");
	need(k->metadata.domain.type == DomainType::WholeStar, "SpinDriveRequiresWholeStar");
	const std::vector<Species> required{{"10", 1, 0}, {"11", 1, 1}, {"0", 0, -1}, {"1", 0, -1}};
	need(k->metadata.species.size() == required.size(), "ForeignStructuralSpecies");
	for (std::size_t i = 0; i < required.size(); ++i)
	{
		const auto &a = k->metadata.species[i], &b = required[i];
		need(a.label == b.label && a.baryon_number == b.baryon_number && a.charge == b.charge, "ForeignStructuralSpecies");
	}
	need(k->metadata.q_normalization == "q=Omega_geom^2; Omega_geom=Omega_phys/c; q in km^-2", "ForeignStructuralQConvention");
	need(k->metadata.profile_units == "r,m:km; epsilon,p:km^-2; n_i=Y_i*n_B:fm^-3; nu:dimensionless", "ForeignStructuralUnits");
	const auto &central_source = k->metadata.sources.front();
	need(std::find(owner->source_files.begin(), owner->source_files.end(), central_source.eos_snapshot.table_path) != owner->source_files.end(), "ForeignStructuralEosAuthority");
	for (const auto &source : k->metadata.sources)
		need(source.table_contents == central_source.table_contents && source.eos_snapshot.identity == central_source.eos_snapshot.identity && source.eos_snapshot.revision == central_source.eos_snapshot.revision && source.eos_snapshot.physical_domain == central_source.eos_snapshot.physical_domain, "MixedStructuralEosAuthority");
	const std::size_t n = z->Channels().size();
	need(vi.size() == n && gn.size() == n && gv.size() == n, "SpinDriveGoalShape");
	RotochemicalSpinDrive out;
	out.chemical_ = z;
	out.structural_ = k;
	out.structural_snapshot_ = structural_text(*k);
	out.structural_sources_snapshot_ = k->metadata.sources;
	out.i_validation_ = vi;
	out.numerical_goal_ = gn;
	out.validation_goal_ = gv;
	const auto physical = k->WholeStarIPhysical(), ke = k->Errors();
	const double ic = AngularVelocity::FromRadPerSecond(1).GeomKmInverse();
	for (auto ch : z->Channels())
	{
		std::string label = ch == ImbalanceChannel::Npe ? "0" : "1";
		auto it = std::find_if(k->metadata.species.begin(), k->metadata.species.end(), [&](const auto &s)
							   { return s.label == label; });
		need(it != k->metadata.species.end(), "StructuralLeptonLabelMissing");
		auto j = it - k->metadata.species.begin();
		out.i_.push_back(physical[j]);
		out.i_error_.push_back(ke[j] * ic * ic);
	}
	const auto zz = z->Values(), ez = z->NumericalError();
	for (std::size_t i = 0; i < n; ++i)
	{
		Sum w, e, v, ar;
		for (std::size_t j = 0; j < n; ++j)
		{
			need(std::isfinite(vi[j]) && vi[j] >= 0, "InvalidStructuralValidationEnvelope");
			w.add(zz[i][j] * out.i_[j]);
			e.add(std::abs(zz[i][j]) * out.i_error_[j] + ez[i][j] * std::abs(out.i_[j]) + ez[i][j] * out.i_error_[j]);
			v.add(std::abs(zz[i][j]) * vi[j]);
			ar.add(32 * u * std::abs(zz[i][j] * out.i_[j]));
		}
		out.w_.push_back(w.value());
		out.arithmetic_.push_back(ar.value());
		out.error_.push_back(e.value() + ar.value());
		out.validation_.push_back(std::nextafter(v.value(), std::numeric_limits<double>::infinity()));
		need(std::isfinite(gn[i]) && gn[i] >= 0 && out.error_.back() <= gn[i], "AccuracyGoalUnmet");
		need(std::isfinite(gv[i]) && gv[i] >= 0 && out.validation_.back() <= gv[i], "StructuralValidationEnvelopeUnmet");
	}
	out.RequireCurrent();
	return out;
}
void RotochemicalSpinDrive::RequireCurrent() const
{
	need(chemical_ && structural_, "MissingSpinDriveDependency");
	chemical_->RequireCurrent();
	coverage(*chemical_->Global()->Lifetime(), *structural_);
	need(structural_->metadata.sources.size() == structural_sources_snapshot_.size(), "StaleStructuralSourceIdentity");
	for (std::size_t i = 0; i < structural_sources_snapshot_.size(); ++i)
	{
		const auto &a = structural_->metadata.sources[i], &b = structural_sources_snapshot_[i];
		need(a.star == b.star && a.profile == b.profile && a.first_order == b.first_order && a.monopole == b.monopole, "StaleStructuralSourceIdentity");
	}
	structural_->RequireCurrent();
	need(structural_text(*structural_) == structural_snapshot_, "StaleStructuralResponse");
}
std::vector<double> RotochemicalSpinDrive::Evaluate(double omega, double dot) const
{
	RequireCurrent();
	need(std::isfinite(omega) && std::isfinite(dot), "NonfiniteSpinInput");
	auto v = w_;
	for (auto &x : v)
	{
		x *= 2 * omega * dot;
		need(std::isfinite(x), "NonfiniteSpinAction");
	}
	return v;
}
} // namespace CompactStar::Analysis
