// Independent background-only characterization, not another accepted G/Z/W object.
void RadialLtildeComparison(const ControlledFixture& f,const RC::FrozenRotochemicalRunContext& reference,const std::filesystem::path& table,const std::filesystem::path& dir) {
 auto star=solve(table.string(),1.10e15,20000,dir);auto rr=star->Profile().GetRadius(),mm=star->Profile().GetMass(),nu=star->Profile().GetMetricNu(),nb=star->Profile().GetBaryonDensity();
 auto value=[&](double x){size_t k=std::upper_bound(rr->Values().begin(),rr->Values().end(),x)-rr->Values().begin();if(k==0)return std::array<double,3>{(*mm)[0]*std::pow(x/(*rr)[0],3),(*nu)[0],(*nb)[0]};if(k==rr->Size())return std::array<double,3>{(*mm)[-1],(*nu)[-1],(*nb)[-1]};double t=(x-(*rr)[k-1])/((*rr)[k]-(*rr)[k-1]);return std::array<double,3>{(*mm)[k-1]+t*((*mm)[k]-(*mm)[k-1]),(*nu)[k-1]+t*((*nu)[k]-(*nu)[k-1]),(*nb)[k-1]+t*((*nb)[k]-(*nb)[k-1])};};
 auto quad=gsl_integration_glfixed_table_alloc(32);
 for(auto process:{RC::UrcaProcess::Me,RC::UrcaProcess::Mmu}) {
  const double onset=process==RC::UrcaProcess::Me?f.provider->NeutronOnsetBaryonDensityFm3():f.provider->MuonOnsetBaryonDensityFm3();double lo=0,hi=(*rr)[-1];for(int j=0;j<80;++j){double mid=(lo+hi)/2;if(value(mid)[2]>onset)lo=mid;else hi=mid;}double edge=(lo+hi)/2;
  std::vector<double> knots{0};for(double x:rr->Values())if(x<edge)knots.push_back(x);knots.push_back(edge);long double total=0;
  for(size_t k=1;k<knots.size();++k)for(int j=0;j<32;++j){double x,w;gsl_integration_glfixed_point(knots[k-1],knots[k],j,&x,&w,quad);auto v=value(x);total+=w*4*M_PI*x*x/std::sqrt(1-2*v[0]/x)*std::exp(-6*v[1]);}
  const double L=double(total)*CompactStar::Units::KM3_TO_CM3*(process==RC::UrcaProcess::Me?1e-51:2e-51),base=reference.Ltilde(process);Near(L,base,5e-3,"predeclared radial20000 Ltilde-only comparison");std::cout<<"RADIAL_LTILDE "<<RC::ProcessIndex(process)<<' '<<L<<" relative "<<std::abs(L-base)/base<<std::endl;
 }
 gsl_integration_glfixed_table_free(quad);
}
