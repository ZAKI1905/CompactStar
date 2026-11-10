#include "process_pool.hpp"
#include <set>
#include <sstream>
using namespace NativeQualification;
void Check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(int argc,char** argv) {
 try {
  Check(argc==2,"fresh output required");std::filesystem::path root=argv[1];
  Check(!std::filesystem::exists(root),"fresh output required");std::filesystem::create_directories(root);
  ProcessPool(73,24,root/"success",[&](std::size_t i,std::size_t slot){
   cpu_set_t a;CPU_ZERO(&a);Check(!sched_getaffinity(0,sizeof a,&a)&&CPU_COUNT(&a)==1,"worker affinity");
   usleep(1000*(7-i%7));std::ofstream f(root/("result-"+std::to_string(i)));f<<i<<' '<<slot<<' '<<sched_getcpu()<<'\n';
  });
  std::set<int> cpus;
  for(std::size_t i=0;i<73;++i){std::ifstream f(root/("result-"+std::to_string(i)));std::size_t got,slot;int cpu;Check(bool(f>>got>>slot>>cpu)&&got==i,"missing result");cpus.insert(cpu);}
  Check(cpus.size()==24,"not all 24 physical workers exercised");
  auto refusal=[&](const char* name,const std::function<void(std::size_t,std::size_t)>& f) {
   bool failed=false;try{ProcessPool(100,2,root/name,f);}catch(const std::exception&){failed=true;}
   Check(failed&&std::filesystem::exists(root/name/"failure.txt"),"worker failure not preserved");
   int status;Check(waitpid(-1,&status,WNOHANG)==-1&&errno==ECHILD,"worker left alive or unreaped");
   std::ifstream in(root/name/"dispatch.tsv");std::string line;std::getline(in,line);std::size_t starts=0;
   while(std::getline(in,line))if(line.find("START")!=std::string::npos)++starts;
   Check(starts==2,"tasks dispatched after early failure");
  };
  refusal("exception",[](std::size_t i,std::size_t){if(i==0)throw std::runtime_error("injected failure");sleep(10);});
  refusal("abrupt-exit",[](std::size_t i,std::size_t){if(i==0)_exit(7);sleep(10);});
  refusal("interrupt",[](std::size_t i,std::size_t){if(i==0)kill(getppid(),SIGTERM);sleep(10);});
  bool missing=false;try{ProcessPool(1,25,root/"invalid",[](auto,auto){});}catch(const std::exception&){missing=true;}
  Check(missing,"unbounded worker count accepted");
  // No numerical calculation: order and reduction are represented by retained task indices.
  std::cout<<"PASS 24 distinct pinned workers, 73 tasks, exception/crash/SIGTERM cleanup, bounded dispatch\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
