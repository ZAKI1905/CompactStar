#include <CompactStar/Analysis/ParticleNumberResponse.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/RelativityUnits.hpp>
#include "../eos/structure1/table.hpp"
#include "particle_number_homogeneous.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <gsl/gsl_integration.h>
using namespace CompactStar;
using namespace CompactStar::Analysis;
using namespace structure1;
const std::vector<Species> species{{"10",1,0},{"11",1,1},{"0",0,-1},{"1",0,-1}};
const NumberDomain whole{DomainType::WholeStar,0,0,"whole star through P=0; canonical finite cut plus explicitly bounded remainder"};
struct Solver : Core::TOVSolver
{
    using Core::TOVSolver::PressureCutoff;
    const auto &labels() const {return eos_tab.extra_labels;}
};
std::shared_ptr<Core::NStar> solve(const std::string &table,double rho,std::size_t resolution,const std::filesystem::path &work,double maximum_radius=70)
{
    require(!std::filesystem::exists(work),"star scratch must be fresh");
    Solver s; s.SetWrkDir(work.string()); s.ImportEOS(table,true); s.SetRadialRes(resolution);
    if(maximum_radius!=70) s.SetMaxRadius(maximum_radius);
    require(rho>s.GetEOSMinEDens()*10 && rho<.999*s.GetEOSMaxEDens(),"central state clamps");
    std::vector<Core::TOVPoint> points;
    require(s.SingleStarSolveToTOVPoints(rho,points)>0 && s.LastSolveStatus()==Core::TOVSolveStatus::SURFACE_REACHED,"Track-R canonical solve failed");
    require(points.back().p==s.PressureCutoff(),"Track-R surface mismatch");
    auto star=std::make_shared<Core::NStar>(points,s.labels());
    require(star->ComputeHartleMonopoleResponse(),"Track-R accepted monopole unavailable");
    return star;
}
NumberTail comparison_tail_bound(const Core::NStar &star);
NumberTail tail_bound(const Core::NStar &star) { return comparison_tail_bound(star); }
NumberInput input(const Core::NStar &star,const std::shared_ptr<NumberEosSource> &source)
{
    std::ostringstream definition; definition<<std::setprecision(17)<<"rho_c=epsilon_phys/c^2 = "
        <<RelativityUnits::EnergyKmMinus2ToMassDensityGcm3((*star.Profile().GetEnergyDensity())[0])
        <<" g/cm^3; canonical achieved central energy in km^-2";
    return {&star,source,species,whole,tail_bound(star),definition.str()};
}
std::shared_ptr<NumberEosSource> source(const std::filesystem::path &table)
{ return std::make_shared<NumberEosSource>(NumberEosSource{"FR2005 human-ratified Structure-1 midpoint free gas","canonical-A1-structure1-table-v1","vacuum,pe,npe,npe-mu; strictly below Sigma ceiling",table.string()}); }
double direct_nodal_count(const Core::NStar &star,int sp,const Core::NStar *reused_metric=nullptr)
{
    const auto &p=star.Profile(); auto R=p.GetRadius(),M=p.GetMass(),NB=p.GetBaryonDensity(),Y=p.GetSpeciesPtr(species[sp].label);
    auto rule=gsl_integration_glfixed_table_alloc(32); long double total=0;
    for(std::size_t i=0;i<R->Size();++i)
    {
        double a=i?(*R)[i-1]:0,b=(*R)[i];
        double n0=i?(*Y)[i-1]*(*NB)[i-1]:(*Y)[0]*(*NB)[0],n1=(*Y)[i]*(*NB)[i];
        for(int j=0;j<32;++j)
        {
            double r,w; gsl_integration_glfixed_point(a,b,j,&r,&w,rule); double t=(r-a)/(b-a);
            double m=i?(*M)[i-1]+t*((*M)[i]-(*M)[i-1]):(*M)[0]*t*t*t;
            if(reused_metric)
            {
                auto rr=reused_metric->Profile().GetRadius(),mm=reused_metric->Profile().GetMass();
                double location=r/(*R)[-1]*(*rr)[-1];
                auto it=std::lower_bound(rr->Values().begin(),rr->Values().end(),location);
                std::size_t k=it-rr->Values().begin();
                if(k==0) m=(*mm)[0]*std::pow(location/(*rr)[0],3);
                else if(k==rr->Size()) m=(*mm)[-1];
                else m=(*mm)[k-1]+(location-(*rr)[k-1])/((*rr)[k]-(*rr)[k-1])*((*mm)[k]-(*mm)[k-1]);
            }
            total+=w*4*M_PI*r*r*(n0+t*(n1-n0))/std::sqrt(1-2*m/r);
        }
    }
    gsl_integration_glfixed_table_free(rule); return double(total)*1e54;
}
// Alternate EOS-knot partition: each original node remains a boundary. Geometry/Hartle
// remain owning-star interpolants; n_i values at inserted knots come from the EOS table.
std::array<double,4> knot_response(const Core::NStar &star,const Table &table,bool refine=true)
{
    const auto &p=star.Profile(); auto R=p.GetRadius(),M=p.GetMass(),NU=p.GetMetricNu(),P=p.GetPressure();
    auto h=star.MonopoleResponse(); const auto &f=star.RotationResponse();
    std::array<double,4> result{};
    for(int sp=0;sp<4;++sp)
    {
        auto Y=p.GetSpeciesPtr(species[sp].label),NB=p.GetBaryonDensity(); std::vector<NumberMeasureSegment> segments;
        std::array<double,3> threshold_terms{}; std::size_t threshold_cells=0;
        auto node=[&](std::size_t i){return NumberMeasureNode{(*R)[i],(*M)[i],(*NU)[i],(*Y)[i]*(*NB)[i],h->xi0_over_Omega2[i],h->m0_over_Omega2[i],f.omega_bar_over_Omega[i]};};
        auto c=node(0); c.r=c.m=c.xi_hat=c.m_hat=0; segments.push_back({c,node(0),false});
        for(std::size_t i=1;i<R->Size();++i)
        {
            auto a=node(i-1),b=node(i); std::vector<NumberMeasureNode> knots{a};
            for(auto it=table.rows.rbegin();refine && it!=table.rows.rend();++it)
            {
                double pressure=RelativityUnits::PressureDynCm2ToKmMinus2(it->pressure_MeV_fm3*Units::MEV_FM3_TO_ERG_CM3);
                if(pressure>=(*P)[i-1] || pressure<=(*P)[i]) continue;
                double t=(pressure-(*P)[i-1])/((*P)[i]-(*P)[i-1]); auto mix=[&](double x,double y){return x+t*(y-x);};
                NumberMeasureNode knot{mix(a.r,b.r),mix(a.m,b.m),mix(a.nu,b.nu),it->number_densities_fm3[sp],mix(a.xi_hat,b.xi_hat),mix(a.m_hat,b.m_hat),mix(a.s,b.s)};
                // A continuous knot mapped to the same representable radius cannot become
                // a fictitious zero-width atom. Retain every distinct interior radius;
                // the original profile endpoints and their exact total Delta n remain.
                if(knot.r>knots.back().r && knot.r<b.r) knots.push_back(knot);
            }
            bool threshold=false;
            for(int onset:{0,3})
            {
                auto y=p.GetSpeciesPtr(species[onset].label);
                threshold=threshold||(((*y)[i-1]>0)!=((*y)[i]>0));
            }
            knots.push_back(b);
            for(std::size_t k=1;k<knots.size();++k)
            {
                segments.push_back({knots[k-1],knots[k],false});
                if(threshold)
                {
                    auto t=IntegrateNumberMeasure({segments.back()},false,false);
                    threshold_terms[0]+=t.density_measure; threshold_terms[1]+=t.metric; threshold_terms[2]+=t.velocity;
                    ++threshold_cells;
                }
            }
        }
        auto v=IntegrateNumberMeasure(segments,true,false); result[sp]=v.density_measure+v.metric+v.velocity+v.boundary;
        auto tail=tail_bound(star);
        std::cout<<"PB6 partition="<<(refine?"EOS-knots":"profile-nodes")<<" nodes="<<R->Size()<<" species="<<sp
                 <<" N="<<v.count<<" A="<<result[sp]<<" density_measure="<<v.density_measure
                 <<" metric="<<v.metric<<" velocity="<<v.velocity<<" outer_boundary="<<v.boundary
                 <<" quadrature_error="<<v.response_error<<" count_quadrature_error="<<v.count_error
                 <<" onset_crossing_cells="<<threshold_cells<<" onset_density_measure="<<threshold_terms[0]
                 <<" onset_metric="<<threshold_terms[1]<<" onset_velocity="<<threshold_terms[2]
                 <<" onset_atoms=0 tail_correction="<<tail.response_correction[sp]
                 <<" comparison_tail_bound="<<tail.response_error[sp]<<" count_tail_bound="<<tail.count_error[sp]
                 <<" valid=true PB13_tail_status=CANDIDATE_PASS\n";
    }
    return result;
}
// Independent finite-q current quadrature: advected density, displaced radial
// geometry, full Jacobian and exact angular Lorentz integral, before expansion.
// It does not invoke the production number measure or define N(q)=N(0)+q*A.
std::array<double,4> finite_current(const Core::NStar &star,double q)
{
    const auto &p=star.Profile(); auto R=p.GetRadius(),M=p.GetMass(),NU=p.GetMetricNu(),NB=p.GetBaryonDensity();
    const auto *h=star.MonopoleResponse(); const auto &f=star.RotationResponse();
    auto rule=gsl_integration_glfixed_table_alloc(16); std::array<long double,4> total{};
    for(std::size_t i=0;i<R->Size();++i)
    {
        double a=i?(*R)[i-1]:0,b=(*R)[i],dr=b-a;
        double xa=i?h->xi0_over_Omega2[i-1]:0,xb=h->xi0_over_Omega2[i];
        double ma=i?(*M)[i-1]:0,mb=(*M)[i];
        for(int j=0;j<16;++j)
        {
            double r,w; gsl_integration_glfixed_point(a,b,j,&r,&w,rule); double t=(r-a)/dr;
            auto mix=[&](double x,double y){return x+t*(y-x);};
            double mass=i?mix(ma,mb):mb*t*t*t,dm=i?(mb-ma)/dr:3*mb*t*t/dr;
            double xi=mix(xa,xb),mh=i?mix(h->m0_over_Omega2[i-1],h->m0_over_Omega2[i]):h->m0_over_Omega2[i]*std::pow(t,5);
            double nu=i?mix((*NU)[i-1],(*NU)[i]):(*NU)[i];
            double nup=i?((*NU)[i]-(*NU)[i-1])/dr:0;
            double s=i?mix(f.omega_bar_over_Omega[i-1],f.omega_bar_over_Omega[i]):f.omega_bar_over_Omega[i];
            double rr=r+q*xi,mm=mass+q*(mh+xi*dm),jacobian=1+q*(xb-xa)/dr;
            double k=q*rr*rr*std::exp(-2*(nu+q*xi*nup))*s*s;
            require(rr>0 && 1-2*mm/rr>0 && jacobian>0 && k>=0 && k<1,"PB11 finite current left perturbative geometry domain");
            double lorentz=k?std::asinh(std::sqrt(k/(1-k)))/std::sqrt(k):1;
            double weight=w*4*M_PI*rr*rr/std::sqrt(1-2*mm/rr)*jacobian*lorentz;
            for(int sp=0;sp<4;++sp)
            {
                auto Y=p.GetSpeciesPtr(species[sp].label);
                double n0=i?(*Y)[i-1]*(*NB)[i-1]:(*Y)[0]*(*NB)[0],n1=(*Y)[i]*(*NB)[i];
                total[sp]+=weight*mix(n0,n1);
            }
        }
    }
    gsl_integration_glfixed_table_free(rule); std::array<double,4> out{};
    for(int i=0;i<4;++i) out[i]=double(total[i])*1e54; return out;
}
void reduction_closure(const std::shared_ptr<NumberEosSource> &src,const std::filesystem::path &dir,
                       const EquilibriumSequenceNumberDerivative &b)
{
    auto central=solve(src->table_path,1.10e15,80000,dir/"central");
    auto start=std::chrono::steady_clock::now(); auto n=ParticleNumbers::Compute(input(*central,src)); n.RequireCurrent();
    std::cout<<"ParticleNumbers_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    start=std::chrono::steady_clock::now(); auto a=FixedCentralEnergyNumberResponse::Compute(input(*central,src)); a.RequireCurrent();
    std::cout<<"A_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    start=std::chrono::steady_clock::now(); auto k=FixedBaryonNumberResponse::Compute(a,b); k.RequireCurrent();
    std::cout<<"K_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    long double AB=static_cast<long double>(a.Values()[0])+a.Values()[1],BB=static_cast<long double>(b.Values()[0])+b.Values()[1];
    std::array<double,4> independent{}; for(int i=0;i<4;++i) independent[i]=double(static_cast<long double>(a.Values()[i])-static_cast<long double>(b.Values()[i])*AB/BB);
    for(int i=0;i<4;++i) require(std::abs(independent[i]-k.Values()[i])<=k.Errors()[i],"PB9 independent raw expansion");
    std::cout<<"PB9 A_B="<<k.A_B<<" B_B="<<k.B_B<<" B_B_error="<<k.B_B_error<<" conditioning="<<k.conditioning
             <<" depsilon_dq="<<k.central_energy_per_q<<" baryon_residual="<<k.baryon_residual<<" budget="<<k.baryon_budget<<'\n';
    std::array<double,4> sign_flip{},wrong_baryon{};
    for(int i=0;i<4;++i) {sign_flip[i]=a.Values()[i]+b.Values()[i]*double(AB/BB); wrong_baryon[i]=a.Values()[i]-b.Values()[i]*a.Values()[1]/b.Values()[1];}
    require(std::abs(sign_flip[0]+sign_flip[1])>k.baryon_budget,"fixed-baryon sign mutation escaped");
    require(std::abs(wrong_baryon[0]+wrong_baryon[1])>k.baryon_budget,"wrong baryon species map escaped");
    std::cout<<"MUTATION FIRED flip fixed-baryon reduction sign\nMUTATION FIRED wrong baryon species map\nPB9 PASS\n";
    // Independent charged-species response uses integration by parts via a small-q
    // current extrapolation; sequence charged counts use the independent nodal rule.
    auto current0=finite_current(*central,0),current1=finite_current(*central,1e-7),current2=finite_current(*central,5e-8);
    std::array<double,4> charged_K{};
    const double h=.00025,eps=b.metadata.central_energy_km_minus2;
    for(int i=0;i<4;++i)
    {
        double ac=2*(current2[i]-current0[i])/5e-8-(current1[i]-current0[i])/1e-7;
        double bc=(8*(direct_nodal_count(*b.contributing_stars[13],i)-direct_nodal_count(*b.contributing_stars[11],i))
                     -(direct_nodal_count(*b.contributing_stars[14],i)-direct_nodal_count(*b.contributing_stars[10],i)))/(12*h*eps);
        charged_K[i]=ac-bc*double(AB/BB);
        require(std::abs(charged_K[i]-k.Values()[i])<=k.Errors()[i],"PB10 independent charged-current/sequence comparison exceeds propagated budget");
    }
    double reconstructed=charged_K[2]+charged_K[3];
    require(std::abs(charged_K[1]-reconstructed)<=k.charge_budget,"PB10 independent neutral reconstruction failed");
    require(std::abs(k.Values()[1]-k.Values()[0]-k.Values()[3])>k.charge_budget,"wrong charge map/order mutation escaped");
    std::cout<<"PB10 raw_charge="<<k.charge_residual<<" budget="<<k.charge_budget<<" independent_charge="<<charged_K[1]-reconstructed
             <<" reconstructed_proton_K="<<reconstructed<<"\nMUTATION FIRED wrong charge map/order\nPB10 PASS\n";
    auto physical=k.WholeStarIPhysical();
    for(int i=0;i<4;++i) std::cout<<"STRUCTURAL species="<<i<<" N="<<n.Values()[i]<<" A="<<a.Values()[i]<<" B="<<b.Values()[i]
        <<" K="<<k.Values()[i]<<" error="<<k.Errors()[i]<<" I_geom="<<k.Values()[i]<<" I_phys="<<physical[i]<<'\n';
    double previous_baryon=0; std::array<double,4> previous_error{};
    start=std::chrono::steady_clock::now();
    for(double q:{1e-6,5e-7,2.5e-7,1.25e-7})
    {
        double rho=1.10e15*(1+k.central_energy_per_q*q/eps);
        auto shifted=solve(src->table_path,rho,80000,dir/("shift-"+std::to_string(int(q*1e10))));
        auto finite=finite_current(*shifted,q); double delta_b=(finite[0]-current0[0])+(finite[1]-current0[1]);
        std::cout<<"PB11 q="<<q<<" Omega_phys="<<std::sqrt(q)*299792.458<<" requested_rho="<<rho
                 <<" achieved_epsilon="<<(*shifted->Profile().GetEnergyDensity())[0]<<" Delta_N_B="<<delta_b
                 <<" baryon_ratio="<<(previous_baryon?delta_b/previous_baryon:0)<<'\n';
        // Order window specified before running: 0.25 for q^2 and 0.5 for q,
        // allowing 20 percent for the independently bounded numerical floor.
        if(previous_baryon) require(delta_b/previous_baryon>.2 && delta_b/previous_baryon<.3,"PB11 baryon closure is not O(q squared)");
        for(int i=0;i<4;++i)
        {
            double quotient=(finite[i]-current0[i])/q,error=quotient-k.Values()[i];
            std::cout<<"PB11 species="<<i<<" quotient="<<quotient<<" quotient_error="<<error
                     <<" error_ratio="<<(previous_error[i]?error/previous_error[i]:0)<<'\n';
            if(previous_error[i]) require(error/previous_error[i]>.4 && error/previous_error[i]<.6,"PB11 species quotient does not approach K linearly");
            previous_error[i]=error;
        }
        previous_baryon=delta_b;
    }
    auto no_shift=finite_current(*central,1e-6);
    require(std::abs(no_shift[0]+no_shift[1]-current0[0]-current0[1])>10*std::abs(previous_baryon),"suppressed PB11 central shift mutation escaped");
    std::cout<<"PB11_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
             <<"\nMUTATION FIRED suppress PB11 central shift\nPB11 PASS\n";
}
void dependency_sweep(const Table &fine_table,const std::shared_ptr<NumberEosSource> &src,const std::filesystem::path &dir)
{
    Solver solver; solver.SetWrkDir((dir/"dEdP-solve").string()); solver.ImportEOS(src->table_path,true); solver.SetRadialRes(80000);
    std::vector<Core::TOVPoint> points;
    require(solver.SingleStarSolveToTOVPoints(1.10e15,points)>0 && solver.LastSolveStatus()==Core::TOVSolveStatus::SURFACE_REACHED,"PB12 canonical input solve");
    Core::NStar original(points,solver.labels()); require(original.ComputeHartleMonopoleResponse(),"PB12 original monopole");
    std::size_t mutations=0;
    for(auto &p:points) if(p.rho>1e-10 && p.rho<1e-6) {p.dedp*=1.00410; ++mutations;}
    require(mutations>0,"PB12 dEdP mutation has no coverage");
    Core::NStar perturbed(points,solver.labels()); require(perturbed.ComputeHartleMonopoleResponse(),"PB12 perturbed monopole");
    const auto &fo=original.RotationResponse(),&fm=perturbed.RotationResponse();
    const auto &ho=*original.MonopoleResponse(),&hm=*perturbed.MonopoleResponse();
    for(std::size_t i=0;i<points.size();++i)
        require(fo.omega_bar_over_Omega[i]==fm.omega_bar_over_Omega[i] && ho.m0_over_Omega2[i]==hm.m0_over_Omega2[i]
                && ho.p0star_over_Omega2[i]==hm.p0star_over_Omega2[i] && ho.xi0_over_Omega2[i]==hm.xi0_over_Omega2[i],"PB12 hidden interior dEdP Hartle dependency");
    auto ao=FixedCentralEnergyNumberResponse::Compute(input(original,src)),am=FixedCentralEnergyNumberResponse::Compute(input(perturbed,src));
    ao.RequireCurrent(); am.RequireCurrent(); require(ao.Values()==am.Values(),"PB12 dEdP-only A changed");
    auto no=ParticleNumbers::Compute(input(original,src)),nm=ParticleNumbers::Compute(input(perturbed,src)); no.RequireCurrent(); nm.RequireCurrent();
    require(no.Values()==nm.Values(),"PB12 dEdP-only count changed");
    std::cout<<"PB12 dEdP_only exact_zero=true mutated_nodes="<<mutations<<" relative_perturbation=0.00410\n";
    TrackRFreeGasThermodynamicProvider eos; std::vector<std::vector<double>> results;
    for(std::size_t resolution:{2048,4096,8192})
    {
        auto table=resolution==8192?fine_table:generate(eos,dir/("eos-"+std::to_string(resolution)+".tsv"),resolution);
        auto s=resolution==8192?src:source(table.file);
        auto star=solve(table.file.string(),1.10e15,80000,dir/("EOS-central-"+std::to_string(resolution)));
        auto a=FixedCentralEnergyNumberResponse::Compute(input(*star,s)); a.RequireCurrent();
        NumberSequenceRecipe recipe{s,species,whole,1.10e15,1.095e15,1.105e15,"smooth midpoint central branch",{.001,.0005,.00025},80000,(dir/("EOS-sequence-"+std::to_string(resolution))).string(),tail_bound};
        recipe.tail_policy_identity="positive-source pe comparison inequalities"; recipe.tail_policy_revision="PB13-comparison-v1";
        auto b=EquilibriumSequenceNumberDerivative::Compute(recipe); b.RequireCurrent();
        auto k=FixedBaryonNumberResponse::Compute(a,b); k.RequireCurrent(); results.push_back(k.Values());
        if(resolution==8192)
        {
            auto km=FixedBaryonNumberResponse::Compute(am,b); km.RequireCurrent();
            require(km.Values()==k.Values(),"PB12 dEdP-only K changed");
        }
        for(int i=0;i<4;++i) std::cout<<"PB12 EOS_resolution="<<resolution<<" species="<<i<<" K="<<k.Values()[i]<<'\n';
    }
    for(std::size_t j=0;j+1<results.size();++j) for(int i=0;i<4;++i)
    {
        double gap=std::abs(results[j][i]/results.back()[i]-1);
        std::cout<<"PB12 EOS_comparison="<<j<<" species="<<i<<" relative="<<gap<<'\n';
        require(gap<=1e-3,"PB12 EOS-value K convergence exceeds 1e-3");
    }
    // Central states on each side of each continuous onset, never on the refusal
    // window itself. These characterize center-source availability independently
    // of the high-density midpoint sequence. The low-density stars have a larger
    // explicitly configured radial domain; no TOV equation or tolerance is changed.
    int index=0;
    for(double onset:{eos.NeutronOnsetBaryonDensityFm3(),eos.MuonOnsetBaryonDensityFm3()})
        for(double factor:{.99,1.01})
        {
            double central_n=onset*factor,rho_c=rho(eos.BarotropeAt(central_n).energy_density_MeV_fm3);
            double maximum_radius=onset==eos.NeutronOnsetBaryonDensityFm3()?50000:70;
            std::cout<<"PB12 center_case="<<index<<" n_c="<<central_n<<" rho_c="<<rho_c<<" max_radius_km="<<maximum_radius<<'\n'<<std::flush;
            auto star=solve(src->table_path,rho_c,80000,dir/("center-"+std::to_string(index++)),maximum_radius);
            auto a=FixedCentralEnergyNumberResponse::Compute(input(*star,src)); a.RequireCurrent();
            auto n=ParticleNumbers::Compute(input(*star,src)); n.RequireCurrent();
            std::cout<<"PB12 center_achieved="<<(*star->Profile().GetEnergyDensity())[0]<<" nodes="<<star->Profile().GetRadius()->Size()
                     <<" central_dEdP="<<(*star->Profile().GetEosDEdP())[0]<<" valid=true onsets_continuous=true atoms=0\n";
        }
    std::cout<<"PB12 PASS\n";
}
NumberTail comparison_tail_bound(const Core::NStar &star)
{
    EnthalpyOracle oracle; const auto &p=star.Profile(); const auto &h=*star.MonopoleResponse(); const auto &f=star.RotationResponse();
    double r=(*p.GetRadius())[-1],m=(*p.GetMass())[-1],ps=(*p.GetPressure())[-1],nu=(*p.GetMetricNu())[-1];
    double H=oracle.H_for_pressure(RelativityUnits::PressureKmMinus2ToDynCm2(ps)); auto endpoint=oracle.at_H(H);
    double e=std::max((*p.GetEnergyDensity())[-1],endpoint.e*oracle.to_geo);
    require(endpoint.ns[0]==0 && endpoint.ns[3]==0,"PB13 comparison tail is only derived for monotone pe exterior");
    double ru=2*m/(1-(1-2*m/r)*std::exp(2*H)),dr=ru-r;
    double mu=m+4*M_PI/3*e*(std::pow(ru,3)-std::pow(r,3)),D=1-2*mu/r,E=std::exp(-2*nu);
    double s=f.omega_bar_over_Omega[-1],sp=f.domega_bar_over_Omega_dr[-1],mh=h.m0_over_Omega2[-1],ph=h.p0star_over_Omega2[-1];
    require(dr>0 && D>0 && s>0 && sp>0 && mh>0 && ph>0,"PB13 positive comparison domain not established");
    double lp=4*M_PI*ru*(e+ps)/D;
    double spu=sp+dr*4*lp/r;
    require(4/ru>lp && s+dr*spu<1,"PB13 first-order comparison envelope not closed");
    double phu=ph+dr*(std::pow(ru,3)*E*spu*spu/12+2*ru*E*(1+ru*spu)/3);
    double gmin=m/(ru*(ru-2*m)),xiu=phu/gmin;
    double mhu=mh+4*M_PI*ru*ru*xiu*e+dr*(std::pow(ru,4)*E*spu*spu/12+8*M_PI/3*std::pow(ru,4)*(e+ps)*E);
    double gmax=(mu+4*M_PI*std::pow(ru,3)*ps)/(r*(r-2*mu));
    double phlower=ph-dr*(mhu*(1+8*M_PI*ru*ru*ps)/std::pow(r-2*mu,2)
        +4*M_PI*(e+ps)*ru*ru*phu/(r-2*mu)+2*ru*ru*E*gmax/3);
    require(phlower>0,"PB13 monopole positive-source comparison not closed");
    double W=4*M_PI*ru*ru/std::sqrt(D);
    NumberTail out{NumberSurface::FiniteCutWithBoundedRemainder,"ADR-0011 preflight section 9 positive-source pe comparison inequalities; canonical A1; PB13 fresh validation",{},{},{},{}};
    for(int i=0;i<4;++i)
    {
        double ns=std::max(endpoint.ns[i],(*p.GetSpeciesPtr(species[i].label))[-1]*(*p.GetBaryonDensity())[-1]);
        double nb=W*ns*dr*1e54,interior=W*xiu*ns*1e54+nb*(mhu/(r-2*mu)+ru*ru*E/3);
        double shell=4*M_PI*r*r/std::sqrt(1-2*m/r)*h.surface_xi0_over_Omega2*(*p.GetSpeciesPtr(species[i].label))[-1]*(*p.GetBaryonDensity())[-1]*1e54;
        out.count_correction.push_back(0);out.count_error.push_back(nb);out.response_correction.push_back(0);out.response_error.push_back(interior+std::abs(shell));
    }
    out.policy_revision="PB13-comparison-v1";
    std::cout<<"PB13 envelope R_upper="<<ru<<" m_upper="<<mu<<" s_prime_upper="<<spu<<" phat_upper="<<phu
             <<" phat_lower="<<phlower<<" xihat_upper="<<xiu<<" mhat_upper="<<mhu<<'\n';
    return out;
}
void surface_validation(const std::shared_ptr<NumberEosSource> &src,const std::filesystem::path &dir,double log_shift_per_q=0)
{
    std::array<double,4> scaleA{},scaleN{}; std::vector<std::array<double,4>> tail_counts;
    std::vector<double> centers{0,-.0005,-.00025,.00025,.0005};
    if(log_shift_per_q) centers={std::log1p(log_shift_per_q*1e-6),std::log1p(log_shift_per_q*5e-7),std::log1p(log_shift_per_q*2.5e-7),std::log1p(log_shift_per_q*1.25e-7)};
    for(std::size_t j=0;j<centers.size();++j)
    {
        double x=centers[j];
        auto star=solve(src->table_path,1.10e15*std::exp(x),80000,dir/("tail-center-"+std::to_string(j+2)));
        auto a=FixedCentralEnergyNumberResponse::Compute(input(*star,src));a.RequireCurrent();
        auto n=ParticleNumbers::Compute(input(*star,src));n.RequireCurrent();
        auto bound=comparison_tail_bound(*star);
        phase5b_oracle::TailContinuation independent;
        auto coarse=independent.solve(*star,1e-9),fine=independent.solve(*star,1e-10);
        tail_counts.push_back(fine.N);
        for(int i=0;i<4;++i)
        {
            if(j==0) {scaleA[i]=std::abs(a.Values()[i]);scaleN[i]=std::abs(n.Values()[i]);}
            double net=fine.A[i]-a.contributions[i][3];
            double nerr=std::abs(fine.N[i]-coarse.N[i]),aerr=std::abs(fine.A[i]-coarse.A[i]);
            std::cout<<"PB13 x="<<x<<" species="<<i<<" R_cut="<<(*star->Profile().GetRadius())[-1]<<" R_vacuum="<<fine.R
                     <<" N_tail="<<fine.N[i]<<" N_bound="<<bound.count_error[i]<<" A_interior_tail="<<fine.A[i]
                     <<" cutoff_atom="<<a.contributions[i][3]<<" net_A_tail="<<net<<" A_bound="<<bound.response_error[i]
                     <<" direct_N_error="<<nerr<<" direct_A_error="<<aerr<<'\n';
            require(fine.N[i]+nerr<=bound.count_error[i] && std::abs(net)+aerr<=bound.response_error[i],"PB13 direct tail outside comparison enclosure");
            require(bound.count_error[i]<=1e-6*scaleN[i] && bound.response_error[i]<=1e-6*scaleA[i],"PB13 comparison bound exceeds 1e-6 scale");
            if(i==0||i==3) require(fine.N[i]==0 && fine.A[i]==0 && bound.count_error[i]==0 && bound.response_error[i]==0,"PB13 extinct species tail nonzero");
        }
    }
    if(!log_shift_per_q) for(int i=0;i<4;++i)
    {
        double coarse=(tail_counts[4][i]-tail_counts[1][i])/.001,fine=(tail_counts[3][i]-tail_counts[2][i])/.0005;
        std::cout<<"PB13 species="<<i<<" dN_tail_dx_coarse="<<coarse<<" dN_tail_dx_fine="<<fine<<" step_error="<<std::abs(coarse-fine)<<'\n';
    }
    std::cout<<"PB13 PASS comparison envelope; factor-16 provisional bound NOT used for acceptance\n";
}
template<class F> void refusal(F f,const char *expected,const char *name)
{
    bool fired=false;
    try { f(); } catch(const std::exception &e) { fired=std::string(e.what()).find(expected)!=std::string::npos; }
    require(fired,name); std::cout<<"MUTATION FIRED "<<name<<'\n';
}
struct MutatedB : EquilibriumSequenceNumberDerivative
{
    explicit MutatedB(const EquilibriumSequenceNumberDerivative &b) {static_cast<EquilibriumSequenceNumberDerivative&>(*this)=b;}
    void zero_baryon() {values_[0]=-values_[1];}
};
void provenance_validation(const std::shared_ptr<NumberEosSource> &src,const std::filesystem::path &dir,
                           EquilibriumSequenceNumberDerivative &b)
{
    auto star=solve(src->table_path,1.10e15,80000,dir/"contract-central");
    auto a=FixedCentralEnergyNumberResponse::Compute(input(*star,src));a.RequireCurrent();
    auto k=FixedBaryonNumberResponse::Compute(a,b);k.RequireCurrent();
    require(b.metadata.sources.size()==15 && b.metadata.tail_inputs.size()==15 && b.metadata.profile_node_counts.size()==15
        && b.metadata.sequence_radial_resolution==80000 && !b.metadata.tail_callback_revision.empty(),"sequence dependency metadata incomplete");
    require(k.metadata.sources.size()==16 && k.metadata.tail_inputs.size()==16,"reduction dependency metadata incomplete");
    double h=.00025,eps=b.metadata.central_energy_km_minus2;
    for(int sp=0;sp<4;++sp)
    {
        double own=(direct_nodal_count(*b.contributing_stars[13],sp)-direct_nodal_count(*b.contributing_stars[11],sp))/(2*h*eps);
        double mutant=(direct_nodal_count(*b.contributing_stars[13],sp,b.contributing_stars[12].get())
            -direct_nodal_count(*b.contributing_stars[11],sp,b.contributing_stars[12].get()))/(2*h*eps);
        require(std::abs(own/b.Values()[sp]-1)<=2e-4,"independent sequence metric control failed");
        require(std::abs(mutant/b.Values()[sp]-1)>2e-4,"central-star metric reuse mutation escaped");
        std::cout<<"MUTATION FIRED reuse central-star metric on neighbor species="<<sp<<" relative="<<std::abs(mutant/b.Values()[sp]-1)<<'\n';
    }
    MutatedB nearzero(b);nearzero.zero_baryon();
    refusal([&]{auto result=FixedBaryonNumberResponse::Compute(a,nearzero);result.RequireCurrent();},"near-zero/ill-conditioned B_B","near-zero B_B");
    auto first=a;first.metadata.sources[0].first_order=&b.contributing_stars[0]->RotationResponse();
    refusal([&]{first.RequireCurrent();},"first-order","stale/mismatched first-order provenance");
    auto monopole=a;monopole.metadata.sources[0].monopole=b.contributing_stars[0]->MonopoleResponse();
    refusal([&]{monopole.RequireCurrent();},"monopole","stale/mismatched monopole provenance");
    auto saved=*src;src->revision+="-mutated";
    refusal([&]{a.RequireCurrent();},"EOS","stale EOS revision");
    refusal([&]{b.RequireCurrent();},"EOS","stale sequence EOS revision");
    *src=saved;
    src->physical_domain+="-mutated";
    refusal([&]{a.RequireCurrent();},"EOS","stale EOS domain");*src=saved;
    std::ifstream file(src->table_path,std::ios::binary);std::string contents((std::istreambuf_iterator<char>(file)),{});file.close();
    {std::ofstream out(src->table_path,std::ios::app);out<<'\n';}
    refusal([&]{b.RequireCurrent();},"EOS","stale exact EOS table bytes");
    {std::ofstream out(src->table_path,std::ios::binary);out<<contents;}
    a.RequireCurrent();b.RequireCurrent();
    // Complete the PB13 tail check at all actually used finite-spin centers before
    // invalidating disposable stars for the final provenance controls.
    surface_validation(src,dir/"shifted-tails",k.central_energy_per_q/eps);
    const_cast<Core::StarProfile&>(b.contributing_stars[0]->Profile()).Touch();
    refusal([&]{b.RequireCurrent();},"stale profile","stale sequence contributing star");
    refusal([&]{k.RequireCurrent();},"stale profile","stale reduced sequence provenance");
    auto preedit_input=input(*star,src);
    const_cast<Core::StarProfile&>(star->Profile()).Touch();
    refusal([&]{a.RequireCurrent();},"stale profile","stale profile provenance");
    auto fresh_attempt=FixedCentralEnergyNumberResponse::Compute(preedit_input);
    require(!fresh_attempt.valid,"recompute with stale owning Hartle inputs accepted");
    std::cout<<"PROVENANCE PASS all sequence contributors and realized tail policies retained; stale inputs refuse\n";
}
void emit_reference(const std::shared_ptr<NumberEosSource> &src,const std::filesystem::path &dir,
                    const EquilibriumSequenceNumberDerivative &b)
{
    auto star=solve(src->table_path,1.10e15,80000,dir/"reference-central");
    auto n=ParticleNumbers::Compute(input(*star,src));n.RequireCurrent();
    auto a=FixedCentralEnergyNumberResponse::Compute(input(*star,src));a.RequireCurrent();
    auto k=FixedBaryonNumberResponse::Compute(a,b);k.RequireCurrent();auto ip=k.WholeStarIPhysical();
    std::ofstream coefficients(dir/"coefficients.tsv");coefficients<<std::setprecision(17)
        <<"species\tlabel\tb\tcharge\tN\tN_error\tA\tA_error\tB\tB_error\tK\tK_error\tI_phys\n";
    const char *names[]={"n","p","e","mu"};
    for(int i=0;i<4;++i) coefficients<<names[i]<<'\t'<<species[i].label<<'\t'<<species[i].baryon_number<<'\t'<<species[i].charge
        <<'\t'<<n.Values()[i]<<'\t'<<n.Errors()[i]<<'\t'<<a.Values()[i]<<'\t'<<a.Errors()[i]
        <<'\t'<<b.Values()[i]<<'\t'<<b.Errors()[i]<<'\t'<<k.Values()[i]<<'\t'<<k.Errors()[i]<<'\t'<<ip[i]<<'\n';
    std::ofstream summary(dir/"summary.tsv");summary<<std::setprecision(17);
    summary<<"central_epsilon_km_minus2\t"<<a.metadata.central_energy_km_minus2
        <<"\nA_B\t"<<k.A_B<<"\nB_B\t"<<k.B_B<<"\nB_B_error\t"<<k.B_B_error<<"\nconditioning\t"<<k.conditioning
        <<"\ndepsilon_dq\t"<<k.central_energy_per_q<<"\ndepsilon_dq_error\t"<<k.central_energy_per_q_error
        <<"\nbaryon_residual\t"<<k.baryon_residual<<"\nbaryon_budget\t"<<k.baryon_budget
        <<"\ncharge_residual\t"<<k.charge_residual<<"\ncharge_budget\t"<<k.charge_budget<<'\n';
    std::ofstream sources(dir/"sources.tsv");sources<<std::setprecision(17)
        <<"role\tprofile_version\tnodes\tepsilon_c\tpressure_cut\tN_tail_error_n\tN_tail_error_p\tN_tail_error_e\tN_tail_error_mu\tA_tail_error_n\tA_tail_error_p\tA_tail_error_e\tA_tail_error_mu\n";
    for(std::size_t j=0;j<k.metadata.sources.size();++j)
    {
        const auto &source=k.metadata.sources[j];const auto &tail=k.metadata.tail_inputs[j];
        sources<<(j?"sequence-"+std::to_string(j-1):"central")<<'\t'<<source.profile_version<<'\t'<<k.metadata.profile_node_counts[j]
            <<'\t'<<(*source.profile->GetEnergyDensity())[0]<<'\t'<<(*source.profile->GetPressure())[-1];
        for(double v:tail.count_error)sources<<'\t'<<v;for(double v:tail.response_error)sources<<'\t'<<v;sources<<'\n';
    }
    std::ofstream stencil(dir/"stencil.tsv");stencil<<std::setprecision(17)<<"step\tspecies\tthree_point\tfive_point\n";
    for(std::size_t j=0;j<b.metadata.steps.size();++j)for(int i=0;i<4;++i)
        stencil<<b.metadata.steps[j]<<'\t'<<names[i]<<'\t'<<b.three_point[j][i]<<'\t'<<b.five_point[j][i]<<'\n';
    require(bool(coefficients)&&bool(summary)&&bool(sources)&&bool(stencil),"candidate reference write failed");
    std::cout<<"CANDIDATE REFERENCE EMITTED; NOT A GOVERNED BASELINE\n";
}
int main(int argc,char **argv)
{
    try
    {
        gsl_set_error_handler_off(); require(argc==3,"usage: phase5b_freegas_validation stage fresh-output-dir");
        std::string stage=argv[1]; auto dir=std::filesystem::absolute(argv[2]); require(!std::filesystem::exists(dir),"validation output must be fresh"); std::filesystem::create_directories(dir);
        std::cout<<std::setprecision(17); TrackRFreeGasThermodynamicProvider eos;
        auto table=generate(eos,dir/"freegas.tsv",8192); auto src=source(table.file);
        if(stage=="PB1")
        {
            // Predeclared cross-method budget: twice the separately measured radial
            // and EOS-table changes, direct-oracle setting change, explicit count-tail
            // enclosure and roundoff. No historical numbers are a target.
            auto coarse_table=generate(eos,dir/"freegas-coarse.tsv",4096);
            auto coarse_source=source(coarse_table.file);
            auto radial=solve(table.file.string(),1.10e15,40000,dir/"radial-coarse");
            auto coarse=solve(coarse_table.file.string(),1.10e15,80000,dir/"EOS-coarse");
            auto fine=solve(table.file.string(),1.10e15,80000,dir/"fine");
            auto nr=ParticleNumbers::Compute(input(*radial,src)); nr.RequireCurrent();
            auto nc=ParticleNumbers::Compute(input(*coarse,coarse_source)); nc.RequireCurrent();
            auto nf=ParticleNumbers::Compute(input(*fine,src)); nf.RequireCurrent();
            phase5b_oracle::Counts oracle;
            auto oc=oracle.solve(1.10e15,1e-10,1e-10),of=oracle.solve(1.10e15,1e-11,1e-11);
            for(int i=0;i<4;++i)
            {
                double same=std::abs(direct_nodal_count(*fine,i)/nf.Values()[i]-1);
                double budget=2*std::abs(nr.Values()[i]-nf.Values()[i])+2*std::abs(nc.Values()[i]-nf.Values()[i])
                    +std::abs(oc[i]-of[i])+nf.Errors()[i]+64*std::numeric_limits<double>::epsilon()*nf.Values()[i];
                double gap=std::abs(nf.Values()[i]-of[i]);
                std::cout<<"PB1 species="<<i<<" N="<<nf.Values()[i]<<" direct_vacuum_N="<<of[i]
                         <<" same_representation_relative="<<same<<" cross_gap="<<gap<<" cross_budget="<<budget<<'\n';
                require(same<=1e-8,"PB1 same-representation independent quadrature failed");
                require(gap<=budget,"PB1 direct-vacuum count exceeds interpolation/oracle budget");
            }
            double nb=direct_nodal_count(*fine,0)+direct_nodal_count(*fine,1);
            double charge=nf.Values()[1]-nf.Values()[2]-nf.Values()[3];
            require(std::abs((nf.Values()[0]+nf.Values()[1])/nb-1)<=1e-8,"PB1 independent baryon count");
            require(std::abs(charge)<=1e-12*(nf.Values()[1]+nf.Values()[2]+nf.Values()[3]),"PB1 raw/neutral-reconstructed charge");
            std::cout<<"PB1 Track-R PASS raw_charge="<<charge<<'\n';
        }
        else if(stage=="PB6")
        {
            std::vector<double> previous; std::vector<std::vector<double>> values;
            for(std::size_t resolution:{10000,20000,40000,80000})
            {
                auto star=solve(table.file.string(),1.10e15,resolution,dir/("radial-"+std::to_string(resolution)));
                auto a=FixedCentralEnergyNumberResponse::Compute(input(*star,src)); a.RequireCurrent();
                double scale=0; for(auto v:a.Values())scale+=std::abs(v);
                std::cout<<"PB6 radial="<<resolution;
                for(auto v:a.Values())std::cout<<'\t'<<v/1e54; std::cout<<'\n';
                auto nodal=knot_response(*star,table,false);
                for(int i=0;i<4;++i) std::cout<<"PB6 diagnostic_sum_roundoff species="<<i<<" difference="<<nodal[i]-a.Values()[i]<<'\n';
                if(!previous.empty())for(int i=0;i<4;++i)require(std::abs(previous[i]-a.Values()[i])<=1e-3*std::abs(a.Values()[i])+1e-8*scale,"PB6 radial partition envelope failed");
                previous=a.Values(); values.push_back(previous);
                if(resolution==80000)
                {
                    auto knots=knot_response(*star,table);
                    for(int i=0;i<4;++i)
                    {
                        double gap=std::abs(knots[i]-previous[i]);
                        std::cout<<"PB6 knot species="<<i<<" relative="<<gap/std::abs(previous[i])<<'\n';
                        require(gap<=1e-3*std::abs(previous[i])+1e-8*scale,"PB6 EOS-knot/threshold partition failed");
                    }
                }
            }
            std::cout<<"PB6 PASS\n";
        }
        else if(stage=="PB12") dependency_sweep(table,src,dir);
        else if(stage=="PB13") surface_validation(src,dir);
        else if(stage=="PB7" || stage=="PB9-11" || stage=="contracts" || stage=="emit")
        {
            NumberSequenceRecipe recipe{src,species,whole,1.10e15,1.095e15,1.105e15,"Structure-1 monotone midpoint smooth central branch",{.001,.0005,.00025},80000,(dir/"sequence").string(),tail_bound};
            auto begin=std::chrono::steady_clock::now(); recipe.tail_policy_identity="positive-source pe comparison inequalities"; recipe.tail_policy_revision="PB13-comparison-v1";
        auto b=EquilibriumSequenceNumberDerivative::Compute(recipe); b.RequireCurrent();
            std::cout<<"sequence_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<'\n';
            for(std::size_t j=0;j<b.five_point.size();++j) for(int i=0;i<4;++i)
                std::cout<<"PB8 h="<<recipe.log_steps[j]<<" species="<<i<<" three_point="<<b.three_point[j][i]<<" five_point="<<b.five_point[j][i]<<'\n';
            for(std::size_t j=0;j<b.metadata.achieved_central_energy_km_minus2.size();++j)
                std::cout<<"PB8 star="<<j<<" achieved_epsilon_c="<<b.metadata.achieved_central_energy_km_minus2[j]
                         <<" profile_version="<<b.metadata.sources[j].profile_version
                         <<" profile_nodes="<<b.contributing_stars[j]->Profile().GetRadius()->Size()<<'\n';
            if(stage=="PB9-11") { reduction_closure(src,dir,b); return 0; }
            if(stage=="contracts") { provenance_validation(src,dir,b); return 0; }
            if(stage=="emit") { emit_reference(src,dir,b); return 0; }
            phase5b_oracle::Homogeneous independent; auto coarse=independent.solve(1.10e15,1e-10,1e-10),fine=independent.solve(1.10e15,1e-11,1e-11);
            for(int i=0;i<4;++i)
            {
                double difference=std::abs(b.Values()[i]/fine.B[i]-1),floor=std::abs(coarse.B[i]/fine.B[i]-1);
                std::cout<<"PB7 species="<<i<<" sequence_B="<<b.Values()[i]<<" homogeneous_B="<<fine.B[i]<<" relative="<<difference<<" oracle_floor="<<floor<<" sequence_error="<<b.Errors()[i]<<'\n';
                require(floor<2e-5,"PB7 homogeneous numerical settings not converged"); require(difference<=2e-4,"PB7 disagreement exceeds 2e-4");
            }
            std::cout<<"PB7 PASS\n";
        }
        else throw std::runtime_error("unknown validation stage");
        return 0;
    }
    catch(const std::exception &e) {std::cerr<<"STOP "<<e.what()<<'\n'; return 1;}
}
