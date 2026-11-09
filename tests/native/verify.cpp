#include "pair.hpp"
int main(int argc,char** argv) {
 try {
    gsl_set_error_handler_off();std::cout<<std::setprecision(17)<<std::unitbuf;
    require(argc==10,"profile thermal entry frozen coefficients archive archive-sha fresh-work output-tsv");
    Inputs inputs;inputs.profile=argv[1];inputs.certificate=inputs.profile/"certificate.txt";
    inputs.thermal=argv[2];inputs.entry=argv[3];inputs.frozen=argv[4];inputs.coefficients=argv[5];inputs.work=argv[8];
    auto sequence=Core::StellarArchive::Load(argv[6],argv[7]);
    MustRefuse([&]{Core::StellarArchive::Load(argv[6],std::string(64,'0'));},"mismatched archive hash");
    MustRefuse([&]{Core::StellarArchive::Load((inputs.work/"missing").string(),argv[7]);},"missing archive");
    auto pair=NativePair(inputs,sequence);
    require(pair.source->OrdinaryContext()==pair.control->OrdinaryContext()&&pair.source->Tangent()==pair.control->Tangent(),"pair backgrounds differ");
    require(pair.control->HasExactZeroChemicalSolution()&&!pair.source->HasExactZeroChemicalSolution(),"source/control isolation contract");
    pair.source->RequireFullCurrent();pair.control->RequireFullCurrent();
    std::ofstream out(argv[9]);out<<std::setprecision(17);
    out<<"role\tt_s\tx\teta_e\teta_mu\txdot\tetadot_e\tetadot_mu\tPdir\tPnet\tCstar\tLgamma\tR18\tBdot\n";
    std::size_t identity=0;
    for(const auto& context:{pair.source,pair.control}) {
      ActualCheckpointContext evaluation(context,++identity);
      require(!evaluation.Breakpoints().empty(),"missing Cstar breakpoint metadata");
      for(double t:{0.0,462269531250.0})for(const auto& eta:std::array<std::array<double,2>,3>{{{{0,0}},{{1e-7,-2e-7}},{{-1e-7,2e-7}}}}) {
        BNV::CheckpointState state{{0,eta[0],eta[1]}},rhs;
        evaluation.Derivative(t,state,rhs);auto d=evaluation.EvaluateDiagnostics(t,state);
        require(d.valid_through_sample&&d.frozen.max_utilization<=1,"frozen validity failed");
        const double G=std::max(1.0,std::abs(d.P_dir_actual_erg_s)+std::abs(d.LH_erg_s)+std::abs(d.DeltaLnu_erg_s)+std::abs(d.Lnu_eq_erg_s)+std::abs(d.Lgamma_erg_s)+std::abs(d.Lother_erg_s)+std::abs(d.L_out_fluid_inf_erg_s));
        require(std::abs(d.R18_residual_erg_s)<=1e-10*G,"R18 energy contract");
        if(context==pair.control) {
          require(d.Bdot_count_s==0&&d.P_dir_actual_erg_s==0&&d.sigma_count_s==std::array<double,2>{{0,0}},"nonzero control source");
          if(eta[0]==0)require(rhs[1]==0&&rhs[2]==0,"zero control equilibrium drift");
        } else require(d.Bdot_count_s==-2.413652026364133e37&&d.P_dir_actual_erg_s>0,"source drive mismatch");
        // Restart state transport at a static state: validates representation, not trajectory accuracy.
        std::stringstream saved;saved<<std::setprecision(17)<<t;for(double v:state)saved<<' '<<v;
        double rt;BNV::CheckpointState restored,rr;saved>>rt;for(auto& v:restored)saved>>v;
        require(rt==t&&restored==state,"restart scalar roundtrip");
        ActualCheckpointContext restart(context,++identity);restart.Derivative(rt,restored,rr);
        require(rr==rhs,"isolated restarted initial RHS mismatch");
        out<<(context==pair.source?"source":"control")<<'\t'<<t;
        for(double v:state)out<<'\t'<<v;for(double v:rhs)out<<'\t'<<v;
        out<<'\t'<<d.P_dir_actual_erg_s<<'\t'<<d.Pnet_erg_s<<'\t'<<d.Cstar_erg_K<<'\t'<<d.Lgamma_erg_s<<'\t'<<d.R18_residual_erg_s<<'\t'<<d.Bdot_count_s<<'\n';
      }
    }
    out.flush();require(bool(out),"static report write failed");
    // Revoking any contributing profile must still invalidate the full pair.
    sequence->contributing_stars.back()->GetSequence(); // mutable accessor version behavior is not assumed
    auto* last=sequence->contributing_stars.back().get();
    last->GetRho_i(std::string("10"));
    MustRefuse([&]{pair.source->RequireFullCurrent();},"mutated neighboring profile");
    std::cout<<"PASS native pair static checks; twelve states; zero integration calls\n";
    return 0;
 } catch(const std::exception& e){std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
