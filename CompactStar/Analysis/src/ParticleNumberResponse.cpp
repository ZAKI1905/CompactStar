#include <CompactStar/Analysis/ParticleNumberResponse.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/Geometry.hpp>
#include <CompactStar/RelativityUnits.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <iomanip>
#include <sstream>
#include <set>
#include <stdexcept>

namespace CompactStar::Analysis
{
namespace
{
constexpr double count_conversion = 1e54;
constexpr double roundoff = 64 * std::numeric_limits<double>::epsilon();
void require(bool b, const std::string &why) { if (!b) throw std::runtime_error(why); }
// Refusal diagnostics only: these helpers never modify or compare integration inputs.
std::uint64_t bits(double x)
{ std::uint64_t out; static_assert(sizeof out==sizeof x); std::memcpy(&out,&x,sizeof x); return out; }
std::string endpoint_text(const char *name,const NumberMeasureNode &p)
{
    std::ostringstream out;
    out<<name<<" r="<<std::setprecision(std::numeric_limits<double>::max_digits10)<<p.r
       <<" density="<<p.density<<" r_hex="<<std::hexfloat<<p.r<<" density_hex="<<p.density
       <<" r_bits=0x"<<std::hex<<bits(p.r)<<" density_bits=0x"<<bits(p.density);
    return out.str();
}
struct ContinuityRefusal : std::runtime_error
{
    std::size_t segment;
    ContinuityRefusal(std::size_t j,const NumberMeasureSegment &previous,const NumberMeasureSegment &current)
        : std::runtime_error(describe(j,previous,current)),segment(j) {}
    static std::string describe(std::size_t j,const NumberMeasureSegment &a,const NumberMeasureSegment &b)
    {
        auto ulp=[](double x,double y) { auto u=bits(x),v=bits(y); return u>v?u-v:v-u; };
        std::ostringstream out;
        out<<"unrepresented gap/jump in number measure; segment_index="<<j
           <<" previous_declared_jump="<<a.declared_jump<<" current_declared_jump="<<b.declared_jump
           <<" boundary="<<(a.declared_jump||b.declared_jump?"declared true jump interface":"ordinary continuous node expected")
           <<"; "<<endpoint_text("previous_right",a.right)<<"; "<<endpoint_text("current_left",b.left);
        if(std::isfinite(a.right.r)&&std::isfinite(b.left.r)&&a.right.r>=0&&b.left.r>=0)
            out<<" r_ULP="<<ulp(a.right.r,b.left.r);
        if(std::isfinite(a.right.density)&&std::isfinite(b.left.density)&&a.right.density>=0&&b.left.density>=0)
            out<<" density_ULP="<<ulp(a.right.density,b.left.density);
        return out.str();
    }
};
struct Sum
{
    double value = 0, correction = 0;
    void add(double x) { double t = value + x; correction += std::abs(value) >= std::abs(x) ? (value-t)+x : (x-t)+value; value=t; }
    double get() const { return value + correction; }
};
double dot(const std::vector<double> &a, const std::vector<Species> &s, bool charge)
{
    require(a.size()==s.size(), "species ordering/size mismatch"); Sum out;
    for (std::size_t i=0;i<a.size();++i) out.add(a[i]*(charge?s[i].charge:s[i].baryon_number));
    return out.get();
}
double error_sum(const std::vector<double> &a, const std::vector<Species> &s, bool charge)
{
    Sum out; for (std::size_t i=0;i<a.size();++i) out.add(a[i]*std::abs(charge?s[i].charge:s[i].baryon_number)); return out.get();
}
std::string contents(const std::string &path)
{
    if (path.empty()) return {};
    std::ifstream in(path, std::ios::binary); require(bool(in), "EOS source file unavailable");
    return std::string(std::istreambuf_iterator<char>(in), {});
}
bool same_eos(const NumberEosSource &a,const NumberEosSource &b)
{ return a.identity==b.identity && a.revision==b.revision && a.physical_domain==b.physical_domain && a.table_path==b.table_path; }
bool same_species(const std::vector<Species> &a,const std::vector<Species> &b)
{
    if (a.size()!=b.size()) return false;
    for(std::size_t i=0;i<a.size();++i) if(a[i].label!=b[i].label || a[i].charge!=b[i].charge || a[i].baryon_number!=b[i].baryon_number) return false;
    return true;
}
NumberMetadata make_metadata(const NumberInput &in, bool response)
{
    require(in.star && in.eos, "missing owning star/EOS provenance");
    require(!in.eos->identity.empty() && !in.eos->revision.empty() && !in.eos->physical_domain.empty(), "missing EOS identity/revision/domain");
    require(!in.central_state_definition.empty() && !in.domain.boundary_definition.empty(), "missing central state/domain semantics");
    require(!in.species.empty(), "missing species map");
    std::set<std::string> labels;
    for(const auto &s:in.species) require(!s.label.empty() && labels.insert(s.label).second && std::isfinite(s.baryon_number) && std::isfinite(s.charge), "invalid/duplicate species map");
    const auto &p=in.star->Profile();
    require(!p.empty() && p.GetRadius() && p.GetPressure() && p.GetEnergyDensity() && p.GetBaryonDensity() && p.GetMass() && p.GetMetricNu(), "incomplete profile");
    NumberMetadata m; m.species=in.species; m.domain=in.domain; m.surface=in.tail.semantics;
    m.surface_authority=in.tail.authority; m.central_state_definition=in.central_state_definition;
    m.tail_inputs.push_back(in.tail); m.profile_node_counts.push_back(p.GetRadius()->Size());
    m.central_energy_km_minus2=(*p.GetEnergyDensity())[0];
    require(m.central_energy_km_minus2>0 && std::isfinite(m.central_energy_km_minus2), "invalid central energy");
    NumberProvenance source; source.star=in.star; source.profile=&p; source.profile_version=p.Version();
    source.eos=in.eos; source.eos_snapshot=*in.eos; source.table_contents=contents(in.eos->table_path);
    if(response || in.tail.semantics==NumberSurface::FiniteCutWithBoundedRemainder)
    {
        source.first_order=&in.star->RotationResponse(); source.monopole=in.star->MonopoleResponse();
        require(source.first_order->MatchesSource(&p,p.Version()) && source.first_order->r_grid==p.GetRadius(), "stale/mismatched first-order Hartle provenance");
        require(source.monopole && source.monopole->MatchesSource(&p,p.Version()) && source.monopole->r_grid==p.GetRadius(), "stale/mismatched monopole Hartle provenance");
        require(source.monopole->I==source.first_order->I, "mismatched Hartle input normalization");
    }
    m.sources.push_back(source); return m;
}
NumberMeasureNode interpolate(const NumberMeasureNode &a,const NumberMeasureNode &b,double t)
{
    auto mix=[&](double x,double y){return x+t*(y-x);};
    NumberMeasureNode v{mix(a.r,b.r),mix(a.m,b.m),mix(a.nu,b.nu),mix(a.density,b.density),mix(a.xi_hat,b.xi_hat),mix(a.m_hat,b.m_hat),mix(a.s,b.s)};
    // Regular-centre extension only; all positive-radius intervals are nodal linear segments.
    if(a.r==0 && a.m==0 && b.r>0) { v.m=b.m*std::pow(t,3); v.m_hat=b.m_hat*std::pow(t,5); }
    return v;
}
std::array<double,4> integrand(const NumberMeasureSegment &s,double t)
{
    auto p=interpolate(s.left,s.right,t); double dr=s.right.r-s.left.r;
    double w=Geometry::ProperVolumeWeight(p.r,p.m);
    double count=w*p.density*dr;
    return {count,-w*p.xi_hat*(s.right.density-s.left.density),
            p.r>0 ? count*p.m_hat/(p.r-2*p.m) : 0,
            count*p.r*p.r*std::exp(-2*p.nu)*p.s*p.s/3};
}
std::array<double,4> quadrature(const NumberMeasureSegment &s,bool fine)
{
    static constexpr double x8[]={.183434642495649805,.525532409916328986,.796666477413626740,.960289856497536232};
    static constexpr double w8[]={.362683783378361983,.313706645877887287,.222381034453374471,.101228536290376259};
    static constexpr double x4[]={.339981043584856265,.861136311594052575};
    static constexpr double w4[]={.652145154862546143,.347854845137453857};
    std::array<Sum,4> sums;
    for(int j=0;j<(fine?4:2);++j) for(double sign:{-1.,1.})
    {
        auto f=integrand(s,(1+sign*(fine?x8[j]:x4[j]))/2);
        for(int k=0;k<4;++k) sums[k].add(f[k]*(fine?w8[j]:w4[j])/2);
    }
    return {sums[0].get(),sums[1].get(),sums[2].get(),sums[3].get()};
}
std::vector<NumberMeasureSegment> segments(const NumberInput &in,std::size_t species,bool response)
{
    const auto &p=in.star->Profile(); auto R=p.GetRadius(), M=p.GetMass(), P=p.GetPressure(), NU=p.GetMetricNu(), NB=p.GetBaryonDensity();
    auto Y=p.GetSpeciesPtr(in.species[species].label);
    const std::size_t n=R->Size(); require(n>=2 && Y && Y->Size()==n && M->Size()==n && P->Size()==n && NU->Size()==n && NB->Size()==n, "missing/misaligned profile species or geometry");
    const auto *h=response?in.star->MonopoleResponse():nullptr;
    const auto *f=response?&in.star->RotationResponse():nullptr;
    if(response) require(h && h->m0_over_Omega2.Size()==n && h->xi0_over_Omega2.Size()==n && f->omega_bar_over_Omega.Size()==n, "misaligned Hartle grid");
    std::vector<NumberMeasureNode> nodes; nodes.reserve(n+1);
    for(std::size_t i=0;i<n;++i)
    {
        require(std::isfinite((*P)[i]) && (*P)[i]>=0 && std::isfinite((*NB)[i]) && (*NB)[i]>=0 && std::isfinite((*Y)[i]) && (*Y)[i]>=0, "invalid pressure/density/fraction");
        if(i) require((*R)[i]>(*R)[i-1] && (*P)[i]<(*P)[i-1], "nonmonotone radial/isobar partition");
        // ADR-0001: reconstruct at NODES before interpolating. Never integrate Y_i.
        nodes.push_back({(*R)[i],(*M)[i],(*NU)[i],(*Y)[i]*(*NB)[i],h?h->xi0_over_Omega2[i]:0,h?h->m0_over_Omega2[i]:0,f?f->omega_bar_over_Omega[i]:0});
    }
    require(in.domain.inner_pressure_km_minus2>=0 && std::isfinite(in.domain.inner_pressure_km_minus2), "invalid inner isobar");
    double outer=nodes.back().r, inner=0;
    auto radius=[&](double pressure)
    {
        require(pressure>(*P)[n-1] && pressure<(*P)[0], "domain boundary crossing/unbracketed isobar");
        std::size_t i=1; while((*P)[i]>pressure) ++i;
        double t=(pressure-(*P)[i-1])/((*P)[i]-(*P)[i-1]); return (*R)[i-1]+t*((*R)[i]-(*R)[i-1]);
    };
    if(in.domain.type==DomainType::ExplicitFixedIsobar)
    {
        require(in.domain.outer_pressure_km_minus2>0 && std::isfinite(in.domain.outer_pressure_km_minus2), "explicit isobar requires positive pressure");
        outer=radius(in.domain.outer_pressure_km_minus2);
    }
    else require(in.domain.outer_pressure_km_minus2==0 && in.domain.inner_pressure_km_minus2==0, "whole-star domain cannot contain an inferred core/shell boundary");
    if(in.domain.inner_pressure_km_minus2>0) inner=radius(in.domain.inner_pressure_km_minus2);
    require(inner<outer, "inverted/empty shell");
    if(nodes.front().r>0)
    {
        auto centre=nodes.front(); centre.r=centre.m=centre.xi_hat=centre.m_hat=0;
        nodes.insert(nodes.begin(),centre);
    }
    // Construct each endpoint once. Original profile nodes are copied verbatim;
    // only explicit domain cuts introduce interpolated endpoints. Adjacent segments
    // then consume the same stored endpoint, without reconstructing t=1 or t=0.
    std::vector<NumberMeasureNode> endpoints;
    for(std::size_t i=1;i<nodes.size();++i)
    {
        const auto &a=nodes[i-1],&b=nodes[i]; if(b.r<=inner || a.r>=outer) continue;
        const double lo=std::max(a.r,inner),hi=std::min(b.r,outer);
        if(endpoints.empty()) endpoints.push_back(lo==a.r?a:interpolate(a,b,(lo-a.r)/(b.r-a.r)));
        endpoints.push_back(hi==b.r?b:interpolate(a,b,(hi-a.r)/(b.r-a.r)));
    }
    std::vector<NumberMeasureSegment> out;
    for(std::size_t i=1;i<endpoints.size();++i) out.push_back({endpoints[i-1],endpoints[i],false});
    require(!out.empty(), "empty integration domain"); return out;
}
NumberMeasureIntegral integrate_species(const NumberInput &in,std::size_t species,bool response)
{
    auto partition=segments(in,species,response);
    try { return IntegrateNumberMeasure(partition,response,response&&in.domain.inner_pressure_km_minus2>0); }
    catch(const ContinuityRefusal &e)
    {
        std::ostringstream out;
        out<<e.what()<<"; species_index="<<species<<" species_label="<<in.species[species].label
           <<" terminal_boundary=false explicit_domain_boundary=false";
        // Whole-star construction has exactly one segment per source-node interval,
        // plus the regular-centre interval when the profile begins at positive r.
        if(in.domain.type==DomainType::WholeStar)
        {
            const auto &p=in.star->Profile();
            auto R=p.GetRadius(); auto Y=p.GetSpeciesPtr(in.species[species].label);
            std::size_t i=e.segment-((*R)[0]>0?1:0);
            NumberMeasureNode canonical; canonical.r=(*R)[i]; canonical.density=(*Y)[i]*(*p.GetBaryonDensity())[i];
            out<<" expected_shared_profile_node_index="<<i<<"; "<<endpoint_text("canonical_source",canonical)
               <<"; construction=shared_canonical_endpoint_list; immutable_refinement_knots=none";
        }
        else out<<" source_node_index=unavailable_for_clipped_partition";
        throw std::runtime_error(out.str());
    }
}
void tail(const NumberInput &in,bool response,std::vector<double> &v,std::vector<double> &err)
{
    if(in.domain.type!=DomainType::WholeStar) return;
    require(!in.tail.authority.empty(), "whole-star surface has no declared authority");
    if(in.tail.semantics==NumberSurface::ExactVacuum)
    {
        require((*in.star->Profile().GetPressure())[-1]==0, "finite table cut is not a P=0 surface"); return;
    }
    const auto &c=response?in.tail.response_correction:in.tail.count_correction;
    const auto &e=response?in.tail.response_error:in.tail.count_error;
    require(!in.tail.policy_revision.empty(), "finite-cut surface lacks tail policy revision");
    require(c.size()==v.size() && e.size()==v.size(), "whole-star tail estimate/bound missing for a species");
    for(std::size_t i=0;i<v.size();++i)
    {
        require(std::isfinite(c[i]) && std::isfinite(e[i]) && e[i]>=0, "invalid tail correction/bound"); v[i]+=c[i]; err[i]+=e[i];
    }
}
void compatible(const NumberResult &a,const NumberResult &b)
{
    a.RequireCurrent(); b.RequireCurrent();
    require(a.metadata.domain==b.metadata.domain && same_species(a.metadata.species,b.metadata.species), "inconsistent domain/species order");
    require(same_eos(a.metadata.sources.front().eos_snapshot,b.metadata.sources.front().eos_snapshot), "inconsistent EOS source/revision/domain");
    require(std::abs(a.metadata.central_energy_km_minus2/b.metadata.central_energy_km_minus2-1)<1e-12, "inconsistent central states");
    require(a.metadata.surface==b.metadata.surface, "inconsistent surface semantics");
}
} // namespace
bool NumberDomain::operator==(const NumberDomain &b) const
{ return type==b.type && outer_pressure_km_minus2==b.outer_pressure_km_minus2 && inner_pressure_km_minus2==b.inner_pressure_km_minus2 && boundary_definition==b.boundary_definition; }
void NumberProvenance::RequireCurrent() const
{
    require(profile && star && &star->Profile()==profile && profile->Version()==profile_version, "stale profile provenance");
    require(eos && same_eos(*eos,eos_snapshot) && contents(eos_snapshot.table_path)==table_contents, "stale EOS revision/domain/content");
    if(first_order) require(&star->RotationResponse()==first_order && first_order->MatchesSource(profile,profile_version), "stale first-order Hartle provenance");
    if(monopole) require(star->MonopoleResponse()==monopole && monopole->MatchesSource(profile,profile_version), "stale monopole Hartle provenance");
}
void NumberResult::RequireCurrent() const
{ require(valid,refusal_reason.empty()?"invalid structural result":refusal_reason); for(const auto &s:metadata.sources) s.RequireCurrent(); }
const std::vector<double> &NumberResult::Values() const { RequireCurrent(); return values_; }
const std::vector<double> &NumberResult::Errors() const { RequireCurrent(); return errors_; }
NumberMeasureIntegral IntegrateNumberMeasure(const std::vector<NumberMeasureSegment> &s,bool outer,bool inner)
{
    require(!s.empty(), "empty number measure"); std::array<Sum,4> totals; Sum ne,ae;
    for(std::size_t j=0;j<s.size();++j)
    {
        const auto &z=s[j];
        for(const auto *p:{&z.left,&z.right})
        {
            require(std::isfinite(p->nu) && std::isfinite(p->density) && p->density>=0 && std::isfinite(p->xi_hat) && std::isfinite(p->m_hat) && std::isfinite(p->s), "invalid number-measure node");
            (void)Geometry::ProperVolumeWeight(p->r,p->m);
        }
        if(j && !(s[j-1].right.r==z.left.r && s[j-1].right.density==z.left.density))
            throw ContinuityRefusal(j,s[j-1],z);
        if(z.declared_jump)
        {
            require(z.left.r==z.right.r && z.left.m==z.right.m && z.left.xi_hat==z.right.xi_hat, "invalid declared jump geometry/displacement");
            totals[1].add(-Geometry::ProperVolumeWeight(z.left.r,z.left.m)*z.left.xi_hat*(z.right.density-z.left.density));
        }
        else
        {
            require(z.right.r>z.left.r, "undeclared zero-length segment/atom");
            auto f=quadrature(z,true),c=quadrature(z,false);
            for(int k=0;k<4;++k) totals[k].add(f[k]);
            ne.add(4*std::abs(f[0]-c[0])+roundoff*std::abs(f[0]));
            for(int k=1;k<4;++k) ae.add(4*std::abs(f[k]-c[k])+roundoff*std::abs(f[k]));
        }
    }
    auto boundary=[](const NumberMeasureNode &p){return Geometry::ProperVolumeWeight(p.r,p.m)*p.density*p.xi_hat;};
    double b=(outer?boundary(s.back().right):0)-(inner?boundary(s.front().left):0);
    // A finite-density terminal atom is this boundary term, once. No vacuum row is appended.
    NumberMeasureIntegral out{totals[0].get()*count_conversion,ne.get()*count_conversion,
        totals[1].get()*count_conversion,totals[2].get()*count_conversion,totals[3].get()*count_conversion,b*count_conversion,
        (ae.get()+roundoff*std::abs(b))*count_conversion};
    return out;
}
ParticleNumbers ParticleNumbers::Compute(const NumberInput &in)
{
    ParticleNumbers out;
    try
    {
        out.metadata=make_metadata(in,false);
        for(std::size_t i=0;i<in.species.size();++i) { auto v=integrate_species(in,i,false); out.values_.push_back(v.count); out.errors_.push_back(v.count_error); }
        tail(in,false,out.values_,out.errors_); out.valid=true; out.RequireCurrent();
    }
    catch(const std::exception &e) { out.valid=false; out.refusal_reason=e.what(); }
    return out;
}
FixedCentralEnergyNumberResponse FixedCentralEnergyNumberResponse::Compute(const NumberInput &in)
{
    FixedCentralEnergyNumberResponse out;
    try
    {
        out.metadata=make_metadata(in,true); out.metadata.differentiation_coordinate="q=Omega_geom^2 at fixed central energy density";
        for(std::size_t i=0;i<in.species.size();++i)
        {
            auto v=integrate_species(in,i,true);
            Sum sum; for(double x:{v.density_measure,v.metric,v.velocity,v.boundary}) sum.add(x);
            out.values_.push_back(sum.get()); out.errors_.push_back(v.response_error);
            out.contributions.push_back({v.density_measure,v.metric,v.velocity,v.boundary});
        }
        tail(in,true,out.values_,out.errors_); out.valid=true; out.RequireCurrent();
    }
    catch(const std::exception &e) { out.valid=false; out.refusal_reason=e.what(); }
    return out;
}
namespace
{
// A private adapter exposes existing canonical metadata, never changes a TOV equation.
class SequenceSolver final : public Core::TOVSolver
{
  public:
    const std::vector<std::string> &labels() const { return eos_tab.extra_labels; }
    double lower() const { return central_eps_floor_factor*eos_tab.eps.front(); }
    double upper() const { return .999*eos_tab.eps.back(); }
    double cutoff() const { return PressureCutoff(); }
};
std::vector<double> derivative_weights(const std::vector<double> &x)
{
    std::vector<double> w(x.size());
    for(std::size_t j=0;j<x.size();++j)
    {
        double denominator=1; for(std::size_t k=0;k<x.size();++k) if(k!=j) denominator*=x[j]-x[k];
        require(std::isfinite(denominator) && denominator!=0, "unbracketed/degenerate achieved central states");
        Sum numerator;
        for(std::size_t k=0;k<x.size();++k) if(k!=j)
        { double term=1; for(std::size_t l=0;l<x.size();++l) if(l!=j && l!=k) term*=-x[l]; numerator.add(term); }
        w[j]=numerator.get()/denominator;
    }
    return w;
}
}
EquilibriumSequenceNumberDerivative EquilibriumSequenceNumberDerivative::Compute(const NumberSequenceRecipe &recipe)
{
    EquilibriumSequenceNumberDerivative out;
    try
    {
        require(recipe.eos && !recipe.eos->table_path.empty(), "complete-star sequence requires an authenticated EOS table");
        require(recipe.tail_for_star && !recipe.branch_definition.empty(), "missing sequence branch/surface policy");
        require(!recipe.tail_policy_identity.empty() && !recipe.tail_policy_revision.empty(), "missing sequence tail callback identity/revision");
        require(recipe.log_steps.size()>=2 && recipe.radial_resolution>0 && !recipe.work_directory.empty(), "missing multiscale sequence controls");
        require(recipe.rho_reference_g_cm3>recipe.minimum_rho_g_cm3 && recipe.rho_reference_g_cm3<recipe.maximum_rho_g_cm3, "reference outside declared sequence branch");
        const auto source_contents=contents(recipe.eos->table_path); const auto source_snapshot=*recipe.eos;
        std::vector<std::vector<double>> errors;
        for(std::size_t level=0;level<recipe.log_steps.size();++level)
        {
            double h=recipe.log_steps[level]; require(std::isfinite(h) && h>0 && (level==0 || h<recipe.log_steps[level-1]), "invalid central-step ladder");
            std::vector<ParticleNumbers> counts; std::vector<double> x;
            for(int k=-2;k<=2;++k)
            {
                double requested=recipe.rho_reference_g_cm3*std::exp(k*h);
                require(requested>recipe.minimum_rho_g_cm3 && requested<recipe.maximum_rho_g_cm3, "central step crosses declared physical branch");
                require(same_eos(*recipe.eos,source_snapshot) && contents(recipe.eos->table_path)==source_contents, "EOS changed during sequence");
                SequenceSolver solver;
                auto work=std::filesystem::path(recipe.work_directory)/("star-"+std::to_string(level)+"-"+std::to_string(k+2));
                require(!std::filesystem::exists(work), "sequence scratch directory is not fresh");
                solver.SetWrkDir(work.string()); solver.ImportEOS(recipe.eos->table_path,true); solver.SetRadialRes(recipe.radial_resolution);
                require(requested>solver.lower() && requested<solver.upper(), "central step would clamp/cross EOS domain");
                std::vector<Core::TOVPoint> points;
                require(solver.SingleStarSolveToTOVPoints(requested,points)>0 && solver.LastSolveStatus()==Core::TOVSolveStatus::SURFACE_REACHED, "neighbor TOV did not reach valid surface");
                require(points.back().p==solver.cutoff(), "neighbor surface does not match canonical cutoff");
                double achieved=points.front().e;
                require(std::isfinite(achieved) && std::abs(std::log(achieved/requested))<h*1e-4, "achieved central state outside differentiation budget");
                x.push_back(std::log(achieved/recipe.rho_reference_g_cm3));
                auto star=std::make_shared<Core::NStar>(points,solver.labels());
                // The tail policy may use the owning accepted Hartle input. Never another star's.
                require(star->ComputeHartleMonopoleResponse(), "neighbor Hartle input unavailable for surface policy");
                NumberInput input{star.get(),recipe.eos,recipe.species,recipe.domain,recipe.tail_for_star(*star),"rho_c=epsilon_phys/c^2 in g/cm^3; canonical fixed-central-energy TOV"};
                auto n=ParticleNumbers::Compute(input); n.RequireCurrent();
                require(n.metadata.tail_inputs.front().policy_revision==recipe.tail_policy_revision,
                        "sequence tail callback revision mismatch");
                if(out.metadata.sources.empty()) { out.metadata=n.metadata; out.metadata.sources.clear(); out.metadata.tail_inputs.clear(); out.metadata.profile_node_counts.clear(); }
                out.metadata.sources.insert(out.metadata.sources.end(),n.metadata.sources.begin(),n.metadata.sources.end());
                out.metadata.tail_inputs.insert(out.metadata.tail_inputs.end(),n.metadata.tail_inputs.begin(),n.metadata.tail_inputs.end());
                out.metadata.profile_node_counts.insert(out.metadata.profile_node_counts.end(),n.metadata.profile_node_counts.begin(),n.metadata.profile_node_counts.end());
                out.metadata.achieved_central_energy_km_minus2.push_back(n.metadata.central_energy_km_minus2);
                out.contributing_stars.push_back(star); counts.push_back(std::move(n));
            }
            require(x[0]<x[1] && x[1]<x[2] && x[2]<x[3] && x[3]<x[4] && x[1]<0 && x[3]>0, "unbracketed/crossed neighboring central states");
            // Differentiate at the central achieved state, using ALL actual abscissae.
            double x0=x[2]; for(auto &v:x) v-=x0;
            auto w5=derivative_weights(x),w3=derivative_weights({x[1],x[2],x[3]});
            double eps=counts[2].metadata.central_energy_km_minus2;
            out.metadata.central_energy_km_minus2=eps;
            std::vector<double> d5,d3,err;
            for(std::size_t i=0;i<recipe.species.size();++i)
            {
                Sum a,b,e; const double centre=counts[2].Values()[i];
                for(int k=0;k<5;++k)
                { a.add(w5[k]*(counts[k].Values()[i]-centre)); e.add(std::abs(w5[k])*(counts[k].Errors()[i]+roundoff*std::abs(counts[k].Values()[i]))); }
                for(int k=0;k<3;++k) b.add(w3[k]*(counts[k+1].Values()[i]-centre));
                d5.push_back(a.get()/eps); d3.push_back(b.get()/eps);
                err.push_back(e.get()/eps+std::abs(d5.back()-d3.back()));
            }
            out.five_point.push_back(d5); out.three_point.push_back(d3); errors.push_back(err);
        }
        const auto last=out.five_point.size()-1;
        out.values_=out.five_point[last]; out.errors_=errors[last];
        for(std::size_t i=0;i<out.values_.size();++i)
        {
            const double change=std::abs(out.values_[i]-out.five_point[last-1][i]);
            const double scale=std::max(std::abs(out.values_[i]),std::abs(out.five_point[last-1][i]));
            require(change<=2e-4*scale+errors[last][i]+errors[last-1][i], "PB8 adjacent derivatives exceed central-step/error budget");
            out.errors_[i]+=change;
            require(scale==0 || out.errors_[i]<=2e-4*scale, "numerically ill-conditioned sequence derivative");
        }
        out.metadata.steps=recipe.log_steps;
        out.metadata.sequence_radial_resolution=recipe.radial_resolution;
        out.metadata.sequence_branch_policy=recipe.branch_definition;
        out.metadata.tail_callback_identity=recipe.tail_policy_identity;
        out.metadata.tail_callback_revision=recipe.tail_policy_revision;
        out.metadata.differentiation_coordinate="x=ln(rho_c/rho_ref); complete independent canonical stars; materialized dN/d epsilon_c [count km^2]";
        out.valid=true; out.RequireCurrent();
    }
    catch(const std::exception &e) { out.valid=false; out.refusal_reason=e.what(); }
    return out;
}
FixedBaryonNumberResponse FixedBaryonNumberResponse::Compute(const FixedCentralEnergyNumberResponse &a,const EquilibriumSequenceNumberDerivative &b)
{
    FixedBaryonNumberResponse out;
    try
    {
        compatible(a,b); require(a.metadata.domain.type==DomainType::WholeStar, "fixed baryon constraint must use the whole star");
        out.metadata=a.metadata; out.metadata.sources.insert(out.metadata.sources.end(),b.metadata.sources.begin(),b.metadata.sources.end());
        out.metadata.tail_inputs.insert(out.metadata.tail_inputs.end(),b.metadata.tail_inputs.begin(),b.metadata.tail_inputs.end());
        out.metadata.profile_node_counts.insert(out.metadata.profile_node_counts.end(),b.metadata.profile_node_counts.begin(),b.metadata.profile_node_counts.end());
        out.metadata.sequence_radial_resolution=b.metadata.sequence_radial_resolution;
        out.metadata.sequence_branch_policy=b.metadata.sequence_branch_policy;
        out.metadata.tail_callback_identity=b.metadata.tail_callback_identity;
        out.metadata.tail_callback_revision=b.metadata.tail_callback_revision;
        out.metadata.steps=b.metadata.steps; out.metadata.achieved_central_energy_km_minus2=b.metadata.achieved_central_energy_km_minus2;
        out.metadata.differentiation_coordinate="q=Omega_geom^2 at fixed whole-star baryon number";
        const auto &s=a.metadata.species;
        out.A_B=dot(a.Values(),s,false); out.B_B=dot(b.Values(),s,false);
        out.B_B_error=error_sum(b.Errors(),s,false);
        std::vector<double> magnitudes; for(double v:b.Values()) magnitudes.push_back(std::abs(v));
        double scale=error_sum(magnitudes,s,false);
        require(scale>0 && std::isfinite(out.B_B) && std::abs(out.B_B)>std::max(100*out.B_B_error,1e-10*scale), "near-zero/ill-conditioned B_B; no fallback division");
        out.conditioning=std::abs(out.B_B)/scale;
        out.central_energy_per_q=-out.A_B/out.B_B;
        double aerr=error_sum(a.Errors(),s,false);
        out.central_energy_per_q_error=(aerr+std::abs(out.central_energy_per_q)*out.B_B_error)/(std::abs(out.B_B)-out.B_B_error);
        for(std::size_t i=0;i<s.size();++i)
        {
            double term=b.Values()[i]*out.central_energy_per_q;
            out.values_.push_back(a.Values()[i]+term);
            out.errors_.push_back(a.Errors()[i]+std::abs(out.central_energy_per_q)*b.Errors()[i]+std::abs(b.Values()[i])*out.central_energy_per_q_error+roundoff*(std::abs(a.Values()[i])+std::abs(term)));
        }
        out.baryon_residual=dot(out.values_,s,false); out.baryon_budget=error_sum(out.errors_,s,false);
        out.charge_residual=dot(out.values_,s,true); out.charge_budget=error_sum(out.errors_,s,true);
        require(std::abs(out.baryon_residual)<=out.baryon_budget, "raw fixed-baryon identity exceeds propagated budget");
        require(std::abs(out.charge_residual)<=out.charge_budget, "raw charge identity exceeds propagated budget");
        out.valid=true; out.RequireCurrent();
    }
    catch(const std::exception &e) { out.valid=false; out.refusal_reason=e.what(); }
    return out;
}
FixedBaryonNumberResponse FixedBaryonNumberResponse::ReduceDomain(const FixedCentralEnergyNumberResponse &a,const EquilibriumSequenceNumberDerivative &b,const FixedBaryonNumberResponse &whole)
{
    FixedBaryonNumberResponse out;
    try
    {
        compatible(a,b); whole.RequireCurrent();
        require(a.metadata.domain.type==DomainType::ExplicitFixedIsobar && whole.metadata.domain.type==DomainType::WholeStar, "domain reduction requires explicit isobar and whole-star constraint");
        require(same_species(a.metadata.species,whole.metadata.species) && same_eos(a.metadata.sources.front().eos_snapshot,whole.metadata.sources.front().eos_snapshot) && std::abs(a.metadata.central_energy_km_minus2/whole.metadata.central_energy_km_minus2-1)<1e-12, "core/whole constraint provenance mismatch");
        out.metadata=a.metadata; out.metadata.sources.insert(out.metadata.sources.end(),b.metadata.sources.begin(),b.metadata.sources.end());
        out.metadata.sources.insert(out.metadata.sources.end(),whole.metadata.sources.begin(),whole.metadata.sources.end());
        for(const auto *metadata:{&b.metadata,&whole.metadata})
        {
            out.metadata.tail_inputs.insert(out.metadata.tail_inputs.end(),metadata->tail_inputs.begin(),metadata->tail_inputs.end());
            out.metadata.profile_node_counts.insert(out.metadata.profile_node_counts.end(),metadata->profile_node_counts.begin(),metadata->profile_node_counts.end());
        }
        out.metadata.steps=b.metadata.steps; out.metadata.achieved_central_energy_km_minus2=b.metadata.achieved_central_energy_km_minus2;
        out.metadata.sequence_radial_resolution=b.metadata.sequence_radial_resolution;
        out.metadata.sequence_branch_policy=b.metadata.sequence_branch_policy;
        out.metadata.tail_callback_identity=b.metadata.tail_callback_identity;
        out.metadata.tail_callback_revision=b.metadata.tail_callback_revision;
        out.central_energy_per_q=whole.central_energy_per_q; out.central_energy_per_q_error=whole.central_energy_per_q_error;
        out.A_B=whole.A_B; out.B_B=whole.B_B; out.B_B_error=whole.B_B_error; out.conditioning=whole.conditioning;
        for(std::size_t i=0;i<a.Values().size();++i)
        {
            double term=b.Values()[i]*out.central_energy_per_q;
            out.values_.push_back(a.Values()[i]+term);
            out.errors_.push_back(a.Errors()[i]+std::abs(out.central_energy_per_q)*b.Errors()[i]+std::abs(b.Values()[i])*out.central_energy_per_q_error+roundoff*(std::abs(a.Values()[i])+std::abs(term)));
        }
        out.charge_residual=dot(out.values_,out.metadata.species,true); out.charge_budget=error_sum(out.errors_,out.metadata.species,true);
        require(std::abs(out.charge_residual)<=out.charge_budget, "domain raw charge identity exceeds propagated budget");
        // Enclosed core/shell baryon response is flux, not a zero-baryon constraint.
        out.baryon_residual=dot(out.values_,out.metadata.species,false); out.baryon_budget=error_sum(out.errors_,out.metadata.species,false);
        out.valid=true; out.RequireCurrent();
    }
    catch(const std::exception &e) { out.valid=false; out.refusal_reason=e.what(); }
    return out;
}
std::vector<double> FixedBaryonNumberResponse::WholeStarIPhysical() const
{
    RequireCurrent(); require(metadata.domain.type==DomainType::WholeStar, "whole-star I mapping requested for a core/shell");
    auto out=values_; const double inverse_c=AngularVelocity::FromRadPerSecond(1).GeomKmInverse();
    for(auto &x:out) x*=inverse_c*inverse_c; return out;
}
std::vector<double> FixedBaryonNumberResponse::FixedIsobarIGeometric(const std::vector<double> &outer_Y,double outer_Q_B,const std::vector<double> &inner_Y,double inner_Q_B) const
{
    RequireCurrent(); require(metadata.domain.type==DomainType::ExplicitFixedIsobar, "core I mapping requires an explicit fixed-isobar domain");
    require(outer_Y.size()==values_.size() && std::isfinite(outer_Q_B), "missing outer-isobar composition/flux");
    const bool shell=metadata.domain.inner_pressure_km_minus2>0;
    require(shell ? inner_Y.size()==values_.size() && std::isfinite(inner_Q_B) : inner_Y.empty() && inner_Q_B==0, "inconsistent inner-isobar flux");
    auto out=values_;
    for(std::size_t i=0;i<out.size();++i)
    {
        require(std::isfinite(outer_Y[i]) && outer_Y[i]>=0 && (!shell || (std::isfinite(inner_Y[i]) && inner_Y[i]>=0)), "invalid boundary composition");
        out[i]-=outer_Y[i]*outer_Q_B; if(shell) out[i]+=inner_Y[i]*inner_Q_B;
    }
    return out;
}
std::vector<double> FixedBaryonNumberResponse::WholeStarEquilibriumNumberRate(double omega,double omega_dot) const
{
    require(std::isfinite(omega) && std::isfinite(omega_dot), "nonfinite physical spin/derivative");
    auto out=WholeStarIPhysical();
    for(auto &v:out) {v*=2*omega*omega_dot; require(std::isfinite(v), "nonfinite equilibrium number driving");}
    return out;
}
} // namespace CompactStar::Analysis
