#pragma once
// Candidate lossless transport of existing producer results. No TOV/Hartle solve on load.
#include <CompactStar/Analysis/ParticleNumberResponse.hpp>
#include <CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp>
#include <fstream>
#include <iomanip>
#include <locale>
#include <type_traits>

namespace CompactStar::Core {
class StellarArchive {
    std::istream* in_ = nullptr;
    std::ostream* out_ = nullptr;
    static void check(bool ok,const char* why) {if(!ok)throw std::runtime_error(why);}
    template<class T> void field(T& x) {
        if constexpr(std::is_enum_v<T>) {int v=static_cast<int>(x);field(v);x=static_cast<T>(v);}
        else {
            if(in_)check(bool(*in_>>x),"truncated stellar archive scalar");else *out_<<x<<'\n';
            if constexpr(std::is_floating_point_v<T>)check(std::isfinite(x),"nonfinite stellar archive");
        }
    }
    void field(std::string& x) {
        if(in_)check(bool(*in_>>std::quoted(x)),"truncated stellar archive string");
        else *out_<<std::quoted(x)<<'\n';
    }
    template<class T> void field(std::vector<T>& x) {
        std::size_t n=x.size();field(n);check(n<=1000000,"stellar archive allocation limit");
        if(in_)x.resize(n);for(auto& v:x)field(v);
    }
    template<class... T> void fields(T&... x) {(field(x),...);}
    void field(Zaki::Vector::DataColumn& c) {
        std::string label=c.Label();auto v=c.Values();fields(label,v);
        if(in_)c=Zaki::Vector::DataColumn(label,v);
    }
    void dataset(Zaki::Vector::DataSet& d) {
        std::size_t n=d.ColCount();field(n);check(n<=64,"stellar column limit");
        if(in_)d.Reserve(n,0);for(std::size_t j=0;j<n;++j)field(d[j]);
    }
    void field(Analysis::Species& x) {fields(x.label,x.baryon_number,x.charge);}
    void field(Analysis::NumberDomain& x) {fields(x.type,x.outer_pressure_km_minus2,x.inner_pressure_km_minus2,x.boundary_definition);}
    void field(Analysis::NumberTail& x) {fields(x.semantics,x.authority,x.count_correction,x.count_error,x.response_correction,x.response_error,x.policy_revision);}
    void field(Analysis::NumberEosSource& x) {fields(x.identity,x.revision,x.physical_domain,x.table_path);}
    void star(NStar& s) {
        auto& p=s.prof_;
        check(p.m_editDepth==0&&!p.m_editDirty,"cannot archive active profile edit");
        dataset(p.radial);dataset(s.B_integrand);
        fields(p.m_version,p.profile_precision,p.seq_point.ec,p.seq_point.m,p.seq_point.r,p.seq_point.pc,p.seq_point.b,p.seq_point.I,
          p.M,p.R,p.z_surf,p.idx_r,p.idx_m,p.idx_nuprime,p.idx_p,p.idx_eps,p.idx_nb,p.idx_nu,p.idx_lambda,
          p.species_labels,p.species_idx,p.eos_dedp_,p.has_eos_dedp_,s.surface_ready,s.MomI);
        auto& r=s.rot_solver;auto& f=r.first_order_response_;auto& m=r.monopole_response_;
        fields(f.source_version,f.I,f.omega_bar_over_Omega,f.domega_bar_over_Omega_dr,f.valid);
        fields(m.source_version,m.m0_over_Omega2,m.p0star_over_Omega2,m.delta_p0_over_Omega2,m.xi0_over_Omega2,
          m.deltaM_over_Omega2,m.surface_shell_mass_over_Omega2,m.surface_xi0_over_Omega2,m.R_surface,m.I,m.valid);
        check(p.size()>2&&p.size()<=80000&&p.has_eos_dedp_&&p.eos_dedp_.Size()==p.size(),"incomplete owned stellar profile");
        check(p.species_labels==std::vector<std::string>({"10","11","0","1"})&&p.species_idx==std::vector<int>({8,9,10,11}),"stellar species layout mismatch");
        check(p.radial.ColCount()==12,"stellar radial schema mismatch");
        for(int k=0;k<8;++k)check(p.GetColumnIndex(static_cast<StarProfile::Column>(k))==k,"stellar column mapping mismatch");
        for(std::size_t k=0;k<p.radial.ColCount();++k)check(p.radial[k].Size()==p.size(),"stellar column length mismatch");
        for(auto* c:{&f.omega_bar_over_Omega,&f.domega_bar_over_Omega_dr,&m.m0_over_Omega2,&m.p0star_over_Omega2,&m.delta_p0_over_Omega2,&m.xi0_over_Omega2})check(c->Size()==p.size(),"Hartle column length mismatch");
        check(f.valid&&m.valid&&f.source_version==p.Version()&&m.source_version==p.Version(),"stale archived Hartle provenance");
        if(in_) {
            // Relocate only addresses. Persisted versions must already match; never re-seal stale data.
            f.source_profile=&p;m.source_profile=&p;
            f.r_grid=static_cast<const StarProfile&>(p).GetRadius();m.r_grid=f.r_grid;
            r.AttachNStar(&s);
        }
        check(s.RotationResponse().MatchesSource(&p,p.Version())&&s.MonopoleResponse()!=nullptr,"restored Hartle owner mismatch");
        auto R=static_cast<const StarProfile&>(p).GetRadius();
        for(std::size_t k=1;k<R->Size();++k)check((*R)[k]>(*R)[k-1],"nonmonotone archived radius");
    }
    void sequence(Analysis::EquilibriumSequenceNumberDerivative& s) {
        std::string schema="CompactStar-owned-sequence-v1-native-binary64";field(schema);
        check(schema=="CompactStar-owned-sequence-v1-native-binary64","stellar archive schema mismatch");
        fields(s.valid,s.refusal_reason,s.values_,s.errors_,s.three_point,s.five_point);
        auto& m=s.metadata;
        fields(m.species,m.domain,m.surface,m.surface_authority,m.central_state_definition,m.tail_inputs,m.profile_node_counts,
          m.sequence_radial_resolution,m.sequence_branch_policy,m.tail_callback_identity,m.tail_callback_revision,
          m.refinement_policy,m.q_normalization,m.profile_units,m.differentiation_coordinate,m.central_energy_km_minus2,
          m.steps,m.achieved_central_energy_km_minus2);
        std::size_t n=s.contributing_stars.size();field(n);check(n==15,"expected exactly fifteen contributors");
        if(in_){s.contributing_stars.resize(n);m.sources.resize(n);}
        check(m.sources.size()==n&&m.tail_inputs.size()==n&&m.profile_node_counts.size()==n&&m.achieved_central_energy_km_minus2.size()==n,"incomplete sequence owners");
        check(m.steps==std::vector<double>({.001,.0005,.00025})&&m.sequence_radial_resolution==80000,"sequence configuration mismatch");
        check(s.values_.size()==4&&s.errors_.size()==4&&s.five_point.size()==3&&s.three_point.size()==3,"sequence derivative shape mismatch");
        for(std::size_t i=0;i<n;++i) {
            if(in_)s.contributing_stars[i]=std::make_shared<NStar>();
            auto& owned=*s.contributing_stars[i];star(owned);auto& p=m.sources[i];
            Analysis::NumberEosSource live;if(!in_){check(bool(p.eos),"missing EOS owner");live=*p.eos;}
            fields(p.profile_version,live,p.eos_snapshot,p.table_contents);
            if(in_){p.star=&owned;p.profile=&owned.Profile();p.first_order=&owned.RotationResponse();p.monopole=owned.MonopoleResponse();p.eos=std::make_shared<const Analysis::NumberEosSource>(live);}
            check(p.profile_version==owned.Profile().Version()&&m.profile_node_counts[i]==owned.Size(),"restored provenance/version mismatch");
        }
        s.RequireCurrent();
    }
  public:
    static void Save(const std::string& path,Analysis::EquilibriumSequenceNumberDerivative& s) {
        s.RequireCurrent();check(!std::filesystem::exists(path),"archive output already exists");
        std::ofstream out(path);out.imbue(std::locale::classic());out<<std::setprecision(17);
        StellarArchive a;a.out_=&out;a.sequence(s);out.flush();check(bool(out),"stellar archive write failed");
    }
    static std::shared_ptr<Analysis::EquilibriumSequenceNumberDerivative> Load(const std::string& path,const std::string& sha256) {
        check(sha256.size()==64&&Physics::Rotochemical::FrozenSource(path).sha256==sha256,"stellar archive SHA-256 mismatch");
        std::ifstream in(path);in.imbue(std::locale::classic());StellarArchive a;a.in_=&in;
        auto s=std::make_shared<Analysis::EquilibriumSequenceNumberDerivative>();a.sequence(*s);
        in>>std::ws;check(in.eof(),"trailing stellar archive data");return s;
    }
};
}
