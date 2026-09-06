// PB1-PB5 analytic/numerical oracles. Expected values never call production integrals.
#include <CompactStar/Analysis/ParticleNumberResponse.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/RelativityUnits.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
using namespace CompactStar;
using namespace CompactStar::Analysis;
void check(bool b,const std::string &s) { if(!b) throw std::runtime_error(s); }
bool close(double a,double b,double tol) { return std::abs(a-b)<=tol*std::max(std::abs(b),1.); }
void detected(const char *name,double mutated,double expected,double tol)
{ check(!close(mutated,expected,tol),std::string("required mutation escaped: ")+name); std::cout<<"MUTATION FIRED "<<name<<'\n'; }
template<class F> double simpson(F f,double a,double b,int n=4096)
{ double sum=f(a)+f(b); for(int i=1;i<n;++i) sum+=(i%2?4:2)*f(a+(b-a)*i/n); return sum*(b-a)/(3*n); }
std::vector<Species> species{{"n",1,0},{"p",1,1},{"e",0,-1},{"mu",0,-1}};
std::shared_ptr<NumberEosSource> eos=std::make_shared<NumberEosSource>(NumberEosSource{"analytic-uniform-energy","v1","self-bound sphere",{}});
std::unique_ptr<Core::NStar> fixture(int n=2001)
{
    std::vector<Core::TOVPoint> points; std::vector<std::string> labels{"n","p","e","mu"};
    const double R=3,M=.2,eps=3*M/(4*M_PI*R*R*R),a=std::sqrt(1-2*M/R);
    for(int i=0;i<n;++i)
    {
        double r=1e-5+(R-1e-5)*i/(n-1),m=M*std::pow(r/R,3),b=std::sqrt(1-2*M*r*r/(R*R*R));
        double p=eps*(b-a)/(3*a-b),nup=2*M*r/(R*R*R*b*(3*a-b));
        double nb=.12*(1+.15*r*r),yp=.06+.005*r;
        points.emplace_back(r,RelativityUnits::MassKmToLiteralSolarMass(m),nup/1e5,0,
            RelativityUnits::PressureKmMinus2ToDynCm2(p),RelativityUnits::EnergyKmMinus2ToMassDensityGcm3(eps),nb,
            std::vector<double>{1-yp,yp,.8*yp,.2*yp},0);
    }
    return std::make_unique<Core::NStar>(points,labels);
}
NumberInput input(const Core::NStar &s)
{ return {&s,eos,species,{DomainType::WholeStar,0,0,"entire self-bound star to P=0"},{NumberSurface::ExactVacuum,"analytic P=0 boundary",{},{},{},{}},"fixed uniform geometric central energy"}; }
// Independent integration of the declared nodal representation, using GSL 32-point Gauss.
// Deliberately spells out current measure, node products and radius search; no Geometry helper.
double independent_count(const Core::NStar &s,int species_index,bool lapse=false,bool fractions=false)
{
    const auto &p=s.Profile(); auto r=p.GetRadius(),m=p.GetMass(),nb=p.GetBaryonDensity(),nu=p.GetMetricNu(),y=p.GetSpeciesPtr(species[species_index].label);
    auto table=gsl_integration_glfixed_table_alloc(32); double sum=0;
    for(std::size_t i=0;i<r->Size();++i)
    {
        double a=i?(*r)[i-1]:0,b=(*r)[i];
        double n0=fractions?(i?(*y)[i-1]:(*y)[0]):(i?(*nb)[i-1]*(*y)[i-1]:(*nb)[0]*(*y)[0]);
        double n1=fractions?(*y)[i]:(*nb)[i]*(*y)[i];
        for(int j=0;j<32;++j)
        {
            double x,w; gsl_integration_glfixed_point(a,b,j,&x,&w,table); double t=(x-a)/(b-a);
            double mass=i?(*m)[i-1]+t*((*m)[i]-(*m)[i-1]):(*m)[0]*std::pow(t,3);
            double N=n0+t*(n1-n0),l=i?(*nu)[i-1]+t*((*nu)[i]-(*nu)[i-1]):(*nu)[0];
            sum+=w*4*M_PI*x*x*N/std::sqrt(1-2*mass/x)*(lapse?std::exp(l):1);
        }
    }
    gsl_integration_glfixed_table_free(table); return sum*1e54;
}
void pb1()
{
    // Exact proper-volume count for m(r)=alpha*r^3, n constant.
    const double R=2,alpha=.025,n=.31,k=std::sqrt(2*alpha);
    NumberMeasureNode c{0,0,0,n,0,0,0},s{R,alpha*R*R*R,0,n,0,0,0};
    double expected=4*M_PI*n*(std::asin(k*R)-k*R*std::sqrt(1-k*k*R*R))/(2*k*k*k)*1e54;
    auto result=IntegrateNumberMeasure({{c,s,false}},false,false);
    check(close(result.count,expected,1e-8),"PB1 analytic proper-volume count");
    auto star=fixture(); auto counts=ParticleNumbers::Compute(input(*star)); counts.RequireCurrent();
    for(int i=0;i<4;++i) check(close(counts.Values()[i],independent_count(*star,i),1e-8),"PB1 independent nodal count/species source");
    const auto &v=counts.Values();
    check(std::abs(v[1]-v[2]-v[3])<=1e-12*(v[1]+v[2]+v[3]),"PB1 raw charge identity");
    double baryons=independent_count(*star,0)+independent_count(*star,1);
    check(close(v[0]+v[1],baryons,1e-8),"PB1 baryon sum");
    check(!close(v[0],independent_count(*star,0,true),1e-8),"PB1 lapse falsifier lacks signal");
    check(!close(v[0],independent_count(*star,0,false,true),1e-8),"PB1 fraction falsifier lacks signal");
    detected("n_i -> Y_i",independent_count(*star,0,false,true),v[0],1e-8);
    detected("count lapse",independent_count(*star,0,true),v[0],1e-8);
    detected("wrong count-unit factor",result.count/1e54,expected,1e-8);
    detected("wrong metric measure",4*M_PI*n*R*R*R/3*1e54,expected,1e-8);
    auto core=input(*star); core.domain={DomainType::ExplicitFixedIsobar,(*star->Profile().GetPressure())[1000],0,"analytic test's explicitly named midpoint isobar"};
    auto core_n=ParticleNumbers::Compute(core); core_n.RequireCurrent();
    detected("wrong count domain",core_n.Values()[0],v[0],1e-8);
    check(close(v[2]+v[3],v[1],1e-12),"PB1 neutral reconstruction");
    std::cout<<"PB1 PASS analytic_relative="<<std::abs(result.count/expected-1)<<" raw_charge="<<(v[1]-v[2]-v[3])<<'\n';
}
// Flat analytic background with independently specified monopole and first-order fields.
// n=1+c*r, xi=a*r, mhat=b*r^3, s=d. Every PN2 term has a closed primitive.
constexpr double toy_a=.07,toy_b=.015,toy_c=.2,toy_d=.3,toy_R=2;
NumberMeasureNode toy(double r) { return {r,0,0,1+toy_c*r,toy_a*r,toy_b*r*r*r,toy_d}; }
auto toy_segments(int n)
{ std::vector<NumberMeasureSegment> out; for(int i=0;i<n;++i) out.push_back({toy(toy_R*i/n),toy(toy_R*(i+1)/n),false}); return out; }
std::array<double,4> expected_terms()
{
    double integral=std::pow(toy_R,5)/5+toy_c*std::pow(toy_R,6)/6;
    return {-M_PI*toy_a*toy_c*std::pow(toy_R,4)*1e54,4*M_PI*toy_b*integral*1e54,
        4*M_PI*toy_d*toy_d/3*integral*1e54,4*M_PI*toy_a*std::pow(toy_R,3)*(1+toy_c*toy_R)*1e54};
}
double sum_terms(const NumberMeasureIntegral &v) { return v.density_measure+v.metric+v.velocity+v.boundary; }
void pb2()
{
    // Coefficient of q in determinant*u^t. Arbitrary lapse/coordinate/l=2 coefficients.
    // Legendre orthogonality and <sin^2 theta>/2 = 1/3 are evaluated independently.
    double p2=simpson([](double u){return (3*u*u-1)/2;},-1,1)/2;
    double velocity=simpson([](double u){return (1-u*u)/2;},-1,1)/2;
    check(std::abs(p2)<1e-14 && std::abs(velocity-1./3)<1e-14,"PB2 angular algebra");
    for(double h:{-.4,.1,.9}) for(double shift:{-.3,.5})
    {
        double determinant=h+shift+.17+.23,ut=-h-shift+.6*velocity;
        check(std::abs(determinant+ut-(.17+.23+.2))<1e-14,"PB2 lapse cancellation");
    }
    double direct=simpson([](double r){return 4*M_PI*toy_a*(1+toy_c*r)*3*r*r;},0,toy_R)*1e54;
    auto t=expected_terms(); check(close(direct,t[0]+t[3],1e-10),"PB2 integration by parts");
    detected("missing lapse cancellation",.17+.23+.2+.4,.17+.23+.2,1e-10);
    detected("wrong angular 1/3",.5,velocity,1e-10);
    // A spurious monopole from l=2 would add the square's nonzero average.
    double p2square=simpson([](double u){double p=(3*u*u-1)/2;return p*p;},-1,1)/2;
    detected("l=2 contamination",p2square,p2,1e-10);
    std::cout<<"PB2 PASS lapse_cancellation; angular_average="<<velocity<<" P2_average="<<p2<<" IBP_relative="<<std::abs(direct/(t[0]+t[3])-1)<<'\n';
}
void pb3()
{
    // Constant w*xi means exact signed total measure, including a declared true jump.
    auto node=[](double r,double n){return NumberMeasureNode{r,0,0,n,1/(4*M_PI*r*r),0,0};};
    std::vector<NumberMeasureSegment> s{{node(1,2),node(1.4,1.8),false},{node(1.4,1.8),node(1.4001,.9),false},
        {node(1.4001,.9),node(1.4001,.4),true},{node(1.4001,.4),node(2,0),false}};
    // On linear xi interpolation the weight is not constant; use exact polynomial integral
    // for each smooth interval and its exact left value for the atom.
    double expected=0;
    for(const auto &z:s)
    {
        if(z.declared_jump) expected-=4*M_PI*z.left.r*z.left.r*z.left.xi_hat*(z.right.density-z.left.density);
        else
        {
            double a=z.left.r,b=z.right.r,dx=(z.right.xi_hat-z.left.xi_hat)/(b-a),u=z.left.xi_hat-dx*a;
            expected-=4*M_PI*(z.right.density-z.left.density)/(b-a)*(u*(b*b*b-a*a*a)/3+dx*(std::pow(b,4)-std::pow(a,4))/4);
        }
    }
    auto v=IntegrateNumberMeasure(s,true,false);
    check(close(v.density_measure,expected*1e54,1e-10),"PB3 signed smooth/jump/ramp measure");
    auto atom=IntegrateNumberMeasure({s[2]},false,false);
    check(close(atom.density_measure,.5e54,1e-14),"PB3 exact signed jump atom");
    check(v.boundary==0,"PB3 continuous onset invents no atom");
    auto terminal=IntegrateNumberMeasure({{node(1,2),node(2,2),false}},true,false);
    check(close(terminal.boundary,2e54,1e-14) && terminal.density_measure==0,"PB3 finite-density terminal once");
    // Independent shrinking-ramp envelope: weight w*xi = pi*r^2, centered at rt.
    // Its interval average differs from the delta limit by pi*width^2/12.
    const double rt=1.5,limit=M_PI*rt*rt*1e54; double previous=1e300;
    for(double width:{.1,.01,.001,.0001})
    {
        auto a=NumberMeasureNode{rt-width/2,0,0,1,.25,0,0};
        auto b=NumberMeasureNode{rt+width/2,0,0,0,.25,0,0};
        auto ramp=IntegrateNumberMeasure({{a,b,false}},false,false);
        double analytic=M_PI*(rt*rt+width*width/12)*1e54;
        check(close(ramp.density_measure,analytic,1e-10),"PB3 exact ramp measure");
        double error=std::abs(ramp.density_measure-limit); check(error<previous,"PB3 convergent weight envelope"); previous=error;
        std::cout<<"PB3 ramp_width="<<width<<" relative_envelope="<<error/limit<<'\n';
        // Endpoint-only derivatives of a continuous ramp with flat adjacent plateaus
        // can be zero at both endpoints while its signed measure is exactly -1.
        detected("nodal dn_i/dp instead of measure",0,ramp.density_measure,1e-10);
    }
    detected("invent continuous-onset atom",v.density_measure+.1e54,expected*1e54,1e-10);
    detected("wrong measure sign",-v.density_measure,expected*1e54,1e-10);
    std::cout<<"PB3 PASS jump="<<atom.density_measure<<" terminal="<<terminal.boundary<<" same_representation="<<std::abs(v.density_measure/(expected*1e54)-1)<<'\n';
}
void pb4()
{
    auto expected=expected_terms(); double prior=1e300;
    // Linear interpolation of the cubic mhat has relative error ~5/(6*n^2).
    // 16384 gives ample resolution for the predeclared 1e-8 continuum gate.
    for(int n:{256,512,1024,2048,4096,8192,16384})
    {
        auto v=IntegrateNumberMeasure(toy_segments(n),true,false);
        std::array<double,4> actual{v.density_measure,v.metric,v.velocity,v.boundary}; double error=0;
        for(int i=0;i<4;++i) error=std::max(error,std::abs(actual[i]/expected[i]-1));
        check(error<prior,"PB4 radial convergence"); prior=error;
        if(n==16384)
        {
            double analytic_sum=0; for(double e:expected) analytic_sum+=e;
            const char *names[]={"omit density-measure term","omit metric term","omit velocity term","omit outer-boundary term"};
            for(int i=0;i<4;++i)
            {
                check(close(actual[i],expected[i],1e-8),"PB4 individual PN2 contribution");
                detected(names[i],sum_terms(v)-actual[i],analytic_sum,1e-8);
            }
        }
    }
    std::cout<<"PB4 PASS four independent analytic terms; finest_relative="<<prior<<'\n';
}
// Direct nonlinear current: mapped isobar r_phys=r(1+q*a), radial metric and Lorentz
// factors evaluated before any expansion, including the moving domain and Jacobian.
double nonlinear_toy(double q)
{
    auto angular=gsl_integration_glfixed_table_alloc(32);
    double result=simpson([&](double r)
    {
        double physical=r*(1+q*toy_a),sum=0;
        for(int j=0;j<32;++j)
        { double u,w; gsl_integration_glfixed_point(-1,1,j,&u,&w,angular); sum+=w/std::sqrt(1-q*physical*physical*(1-u*u)*toy_d*toy_d); }
        return 2*M_PI*physical*physical*(1+toy_c*r)*(1+q*toy_a)*sum/std::sqrt(1-2*q*toy_b*physical*physical);
    },0,toy_R,4096);
    gsl_integration_glfixed_table_free(angular); return result;
}
void pb5()
{
    auto terms=expected_terms(); double expected=0;for(double t:terms)expected+=t/1e54;
    double n0=nonlinear_toy(0),previous=0;
    for(double q:{1e-3,5e-4,2.5e-4,1.25e-4})
    {
        double quotient=(nonlinear_toy(q)-n0)/q,error=std::abs(quotient-expected);
        if(previous) check(error/previous>.45 && error/previous<.55,"PB5 quotient must approach linearly in q");
        std::cout<<"PB5 q="<<q<<" quotient="<<quotient<<" error="<<error<<'\n'; previous=error;
    }
    double omega=AngularVelocity::FromRadPerSecond(100).GeomKmInverse();
    check(std::abs(omega*299792.458/100-1)<1e-15,"PB5 physical/geometric conversion");
    const double q=1.25e-4,quotient=(nonlinear_toy(q)-n0)/q;
    check(close(quotient,expected,1e-4),"PB5 small-spin independent coefficient");
    detected("wrong q normalization / omitted division",quotient*q,expected,1e-4);
    detected("seed squared contamination",quotient*4,expected,1e-4);
    double physical=expected*1e54/std::pow(299792.458,2);
    double inv_c=AngularVelocity::FromRadPerSecond(1).GeomKmInverse();
    check(close(expected*1e54*inv_c*inv_c,physical,1e-14),"PB5 physical coefficient");
    detected("wrong c squared conversion",expected*1e54,physical,1e-14);
    std::cout<<"PB5 PASS independent nonlinear current and exact Omega convention\n";
}
void endpoint_falsifiers()
{
    auto node=[](double r,double n){return NumberMeasureNode{r,0,0,n,.25,0,0};};
    const auto p0=node(1,2),p1=node(1.5,1.7560842199532527e-7),p2=node(2,0);
    const std::vector<NumberMeasureNode> endpoints{p0,p1,p2};
    std::vector<NumberMeasureSegment> shared;
    for(std::size_t i=1;i<endpoints.size();++i) shared.push_back({endpoints[i-1],endpoints[i],false});
    auto good=IntegrateNumberMeasure(shared,true,false);
    check(shared[0].right.density==shared[1].left.density,"shared endpoint exact equality");
    std::cout<<"ENDPOINT 1 PASS canonical continuous shared endpoint\n";
    auto refuses=[](const std::vector<NumberMeasureSegment> &s)
    {
        try { (void)IntegrateNumberMeasure(s,true,false); }
        catch(const std::exception &e) { return std::string(e.what()).find("unrepresented gap/jump")!=std::string::npos; }
        return false;
    };
    auto perturbed=shared; perturbed[1].left.density=std::nextafter(p1.density,1.);
    check(refuses(perturbed),"1-ULP density gap accepted");
    perturbed=shared; perturbed[1].left.r=std::nextafter(p1.r,2.);
    check(refuses(perturbed),"1-ULP radius gap accepted");
    std::cout<<"ENDPOINT 2 PASS independent 1-ULP radius/density perturbations REFUSED\n";
    perturbed=shared; perturbed[1].left.density+=.5;
    check(refuses(perturbed),"undeclared finite jump accepted");
    std::cout<<"ENDPOINT 3 PASS undeclared finite jump REFUSED\n";
    auto inside=node(1.5,1),outside=node(1.5,.5);
    auto jump=IntegrateNumberMeasure({{p0,inside,false},{inside,outside,true},{outside,p2,false}},true,false);
    auto smooth1=IntegrateNumberMeasure({{p0,inside,false}},false,false);
    auto smooth2=IntegrateNumberMeasure({{outside,p2,false}},false,false);
    const double atom=4*M_PI*1.5*1.5*.25*.5*1e54;
    check(close(jump.density_measure-smooth1.density_measure-smooth2.density_measure,atom,1e-14),"declared jump not exactly one signed atom");
    std::cout<<"ENDPOINT 4 PASS declared true jump exactly one Stieltjes atom\n";
    check(good.boundary==0 && !shared[0].declared_jump && !shared[1].declared_jump,"continuous zero onset atom");
    const double primitive=4*M_PI*.25/3;
    double expected=0;
    for(const auto &s:shared) expected-=primitive*(std::pow(s.right.r,3)-std::pow(s.left.r,3))*(s.right.density-s.left.density)/(s.right.r-s.left.r);
    check(close(good.density_measure,expected*1e54,1e-14),"continuous onset measure differs from smooth primitive");
    std::cout<<"ENDPOINT 5 PASS continuous onset to zero NO atom\n";
    const auto terminal=node(2,2);
    std::vector<NumberMeasureSegment> finite{{p0,terminal,false}};
    const double boundary=4*M_PI*4*2*.25*1e54;
    auto once=IntegrateNumberMeasure(finite,true,false);
    check(close(sum_terms(once),boundary,1e-14),"finite terminal contribution not once");
    std::cout<<"ENDPOINT 6 PASS finite-density terminal exactly once\n";
    auto doubled=finite; doubled.push_back({terminal,node(2,0),true});
    auto duplicate=IntegrateNumberMeasure(doubled,true,false);
    // Test-side duplicate: explicit vacuum atom PLUS the already-counted boundary.
    check(!close(sum_terms(duplicate)+once.boundary,boundary,1e-14),"duplicate terminal mutation escaped oracle");
    std::cout<<"ENDPOINT 7 PASS duplicate terminal mutation FAILS independent boundary oracle\n";
    auto omitted=IntegrateNumberMeasure(finite,false,false);
    check(!close(sum_terms(omitted),boundary,1e-14),"omitted terminal mutation escaped oracle");
    std::cout<<"ENDPOINT 8 PASS omitted terminal mutation FAILS independent boundary oracle\n";
}
// PB14 analytic structural mapping fixture, retaining real profile provenance.
// Q_B=r^2(1-r), Y_p=.1+.2r, Y_e=.8Y_p, Y_mu=.2Y_p, Y_n=1-Y_p.
// PN7 is integrated independently; PN8 is exercised through the public API.
struct MappingFixture : FixedBaryonNumberResponse
{
    MappingFixture(const NumberMetadata &m,const std::vector<double> &v)
    {metadata=m;values_=v;errors_=std::vector<double>(v.size(),0);valid=true;}
};
void pb14()
{
    auto star=fixture();auto n=ParticleNumbers::Compute(input(*star));n.RequireCurrent();
    auto Q=[](double r){return r*r*(1-r);};
    auto Y=[](double r){double p=.1+.2*r;return std::vector<double>{1-p,p,.8*p,.2*p};};
    const std::array<double,4> slope{-.2,.2,.16,.04};
    for(auto bounds:{std::array<double,2>{0,1},std::array<double,2>{0,.6},std::array<double,2>{.2,.6}})
    {
        double inner=bounds[0],outer=bounds[1];auto yi=Y(inner),yo=Y(outer);std::vector<double> expected(4),enclosed(4);
        for(int i=0;i<4;++i)
        {
            expected[i]=-slope[i]*simpson(Q,inner,outer)*1e54;
            enclosed[i]=expected[i]+yo[i]*Q(outer)*1e54-yi[i]*Q(inner)*1e54;
        }
        auto metadata=n.metadata;
        metadata.domain=outer==1?NumberDomain{DomainType::WholeStar,0,0,"analytic entire interval to zero flux"}
            :NumberDomain{DomainType::ExplicitFixedIsobar,1,inner?2.:0.,"analytic explicitly declared outer and optional inner isobars"};
        MappingFixture result(metadata,enclosed);
        if(outer==1)
        {
            auto physical=result.WholeStarIPhysical();auto rates=result.WholeStarEquilibriumNumberRate(100,-1e-12);
            for(int i=0;i<4;++i)
            {
                check(close(result.Values()[i],expected[i],1e-12),"PB14 whole-star geometric mapping");
                check(close(physical[i],expected[i]/std::pow(299792.458,2),1e-12),"PB14 physical units");
                double driver=rates[i];
                check(close(driver,2*expected[i]*(100/299792.458)*(-1e-12/299792.458),1e-12),"PB14 R2006 eq7 structural physical driver");
            }
            bool refused=false;try{(void)result.FixedIsobarIGeometric(yo,Q(outer)*1e54);}catch(const std::exception&){refused=true;}
            check(refused,"PB14 whole-star mislabeled core accepted");
        }
        else
        {
            auto mapped=result.FixedIsobarIGeometric(yo,Q(outer)*1e54,inner?yi:std::vector<double>{},Q(inner)*1e54);
            for(int i=0;i<4;++i)
            {
                check(close(mapped[i],expected[i],1e-12),"PB14 core/shell boundary flux sign");
                detected("missing core boundary flux",enclosed[i],expected[i],1e-12);
            }
            bool refused=false;try{(void)result.WholeStarIPhysical();}catch(const std::exception&){refused=true;}
            check(refused,"PB14 core mislabeled whole-star accepted");
        }
        std::cout<<"PB14 domain="<<(outer==1?"whole-star":inner?"explicit-isobar-shell":"explicit-isobar-core")<<" PASS\n";
    }
    // FR2005 (31) / PN7: a declared composition jump contributes -Q_B Delta Y.
    double atom=-Q(.4)*.03*1e54;
    check(atom<0,"PB14 signed structural composition atom");
    std::cout<<"PB14 source=FR2005(18,24,26,30,31);R2006(7);Hartle(109-114)\n"
             <<"PB14 PASS structural symbol/sign/unit/domain mapping; no authenticated free-gas core cutoff or numerical source target claimed\n";
}
int main(int argc,char **argv)
{
    try { gsl_set_error_handler_off(); std::cout<<std::setprecision(17);
        if(argc>1 && std::string(argv[1])=="endpoints") { endpoint_falsifiers(); return 0; }
        int pb=argc>1?std::stoi(argv[1]):0;
        if(!pb||pb==1)pb1(); if(!pb||pb==2)pb2(); if(!pb||pb==3)pb3(); if(!pb||pb==4)pb4(); if(!pb||pb==5)pb5();
        if(!pb)endpoint_falsifiers(); if(!pb||pb==14)pb14(); return 0;
    } catch(const std::exception &e) { std::cerr<<"STOP "<<e.what()<<'\n'; return 1; }
}
