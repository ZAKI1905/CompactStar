// Rehydrate a completed accepted history; never claims adaptive process restart.
// Inputs are authenticated by continuation.py before invoking this executable.
BNV::MainTrajectoryResult ResumeMain(Ctx p,double a,State initial,double end,Tol tol,double cap,
 const std::vector<double>& times,const std::filesystem::path& dir) {
 if(!std::filesystem::exists(dir/"main.accepted_states.tsv"))return native_main(p,a,initial,end,tol,cap,times,dir);
 require(a==0&&initial==State{{0,0,0}},"only original-t0 retained histories supported");
 BNV::MainTrajectoryResult main;main.start_time_s=a;main.initial_state=initial;main.final_time_s=end;
 main.configuration.relative=tol.r;main.configuration.absolute=tol.a;main.configuration.initial_step_s=1;main.configuration.maximum_step_s=cap;
 main.provenance={p->RunCardIdentity(),"native snapshot mathematical source","full/cheap unchanged guards","RF-native-passive"};
 std::ifstream states(dir/"main.accepted_states.tsv"),events(dir/"events.tsv"),internal(dir/"main.internal_steps.tsv");
 require(bool(states)&&bool(events)&&bool(internal),"retained history missing");
 std::string line,event,detail;std::getline(states,line);std::getline(events,event);std::getline(internal,detail);
 State left=initial;double previous=a;
 while(std::getline(states,line)) {
  BNV::AcceptedStepRecord s;std::istringstream r(line);
  r>>s.ordinal>>s.t_left_s>>s.t_right_s>>s.accepted_step_s>>s.suggested_next_h_s>>s.cumulative_rejected;
  for(auto& v:s.y_right)r>>v;s.y_left=left;left=s.y_right;
  require(bool(r)&&s.ordinal==main.accepted_steps.size()+1&&s.t_left_s==previous&&s.t_right_s>s.t_left_s,"bad retained history");previous=s.t_right_s;
  require(bool(std::getline(events,event))&&bool(std::getline(internal,detail)),"missing retained metadata");
  std::size_t ordinal;std::istringstream e(event);e>>ordinal>>s.event.time>>s.event.knot>>s.event.direction>>s.event.residual>>s.event.time_low>>s.event.time_high;
  for(auto& v:s.local_rk_error)e>>v;require(bool(e)&&ordinal==s.ordinal,"bad retained event");
  if(s.event.knot>=0){s.event.state=s.y_right;s.event.rk_error=s.local_rk_error;s.event.time_ulp=std::nextafter(s.event.time,INFINITY)-s.event.time;}
  s.source_identity=main.provenance.source_identity;s.currentness_identity=main.provenance.currentness_identity;
  std::istringstream d(detail);double ceiling,t0,t1,h,next,x0,x1;std::size_t cp,rejected;
  d>>ordinal>>cp>>ceiling>>t0>>t1>>h>>next>>rejected>>x0>>x1>>s.cstar_cell_left>>s.cstar_cell_right>>s.cstar_knots_crossed;
  require(bool(d)&&ordinal==s.ordinal&&t0==s.t_left_s&&t1==s.t_right_s&&next==s.suggested_next_h_s&&rejected==s.cumulative_rejected&&x0==s.y_left[0]&&x1==s.y_right[0],"retained internal mismatch");
  main.accepted_steps.push_back(s);
 }
 require(!main.accepted_steps.empty()&&previous==end,"incomplete retained main");
 require(!std::getline(events,event)&&!std::getline(internal,detail),"extra retained metadata");
 main.final_state=main.accepted_steps.back().y_right;
 main.statistics.accepted_steps=main.accepted_steps.size();main.statistics.rejected_steps=main.accepted_steps.back().cumulative_rejected;
 require(main.statistics.rejected_steps<=.1*(main.statistics.accepted_steps+main.statistics.rejected_steps),"retained rejection guard failed");
 std::cout<<"REUSED_ACCEPTED_HISTORY "<<dir<<" steps="<<main.accepted_steps.size()<<'\n';return main;
}
