#pragma once
// Linux qualification harness only. One forked, pinned process per worker;
// immutable numerical inputs are inherited by copy-on-write, never shared caches.
#ifdef __linux__
#include <sys/socket.h>
#include <sys/wait.h>
#include <sys/prctl.h>
#include <sched.h>
#include <poll.h>
#include <unistd.h>
#include <signal.h>
#include <cerrno>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>

namespace NativeQualification {
inline volatile sig_atomic_t pool_signal=0;
inline void PoolSignal(int value){pool_signal=value;}
inline std::vector<int> PhysicalCpus() {
 cpu_set_t allowed;CPU_ZERO(&allowed);
 if(sched_getaffinity(0,sizeof allowed,&allowed))throw std::runtime_error("cannot read execution affinity");
 std::map<std::pair<int,int>,int> unique;
 for(int cpu=0;cpu<CPU_SETSIZE;++cpu)if(CPU_ISSET(cpu,&allowed)) {
  const std::string p="/sys/devices/system/cpu/cpu"+std::to_string(cpu)+"/topology/";
  int package,core;std::ifstream a(p+"physical_package_id"),b(p+"core_id");
  if(!(a>>package)||!(b>>core))throw std::runtime_error("missing physical CPU topology");
  unique.emplace(std::make_pair(package,core),cpu);
 }
 std::vector<int> result;for(const auto& item:unique)result.push_back(item.second);return result;
}
struct PoolReply {std::uint64_t index;std::uint64_t ok;};
inline void Send(int fd,const void* p,std::size_t n) {
 ssize_t sent;do{sent=send(fd,p,n,MSG_NOSIGNAL);}while(sent<0&&errno==EINTR&&!pool_signal);
 if(sent!=ssize_t(n))throw std::runtime_error("worker IPC send failed");
}
inline void Receive(int fd,void* p,std::size_t n) {
 ssize_t got;do{got=recv(fd,p,n,0);}while(got<0&&errno==EINTR&&!pool_signal);
 if(got!=ssize_t(n))throw std::runtime_error("worker exited or IPC truncated");
}
struct PoolChild {pid_t pid=-1;int fd=-1;std::uint64_t task=0;bool busy=false;};
inline void StopChildren(std::vector<PoolChild>& children) noexcept {
 for(auto& c:children){if(c.fd>=0){close(c.fd);c.fd=-1;}if(c.pid>0)kill(c.pid,SIGTERM);}
 for(int trial=0;trial<50;++trial) {
  bool pending=false;
  for(auto& c:children)if(c.pid>0){int status;const auto r=waitpid(c.pid,&status,WNOHANG);if(r==c.pid||(r<0&&errno==ECHILD))c.pid=-1;else pending=true;}
  if(!pending)return;usleep(10000);
 }
 for(auto& c:children)if(c.pid>0){kill(c.pid,SIGKILL);while(waitpid(c.pid,nullptr,0)<0&&errno==EINTR){}c.pid=-1;}
}
// callback writes its complete per-task evidence before reporting success/failure.
// Aggregate/reduction is exclusively the caller's ordered serial operation.
inline void ProcessPool(std::size_t tasks,std::size_t workers,const std::filesystem::path& evidence,
                        const std::function<void(std::size_t,std::size_t)>& callback) {
 if(!tasks||!workers||workers>24)throw std::runtime_error("invalid bounded worker request");
 const auto cpus=PhysicalCpus();if(cpus.size()<workers)throw std::runtime_error("insufficient distinct allocated physical cores");
 std::filesystem::create_directories(evidence);
 if(std::filesystem::exists(evidence/"dispatch.tsv"))throw std::runtime_error("process pool output already exists");
 std::ofstream dispatch(evidence/"dispatch.tsv");dispatch<<"index\tworker\tpid\tcpu\tevent\n";dispatch.flush();
 struct sigaction action{},old_term{},old_int{};action.sa_handler=PoolSignal;sigemptyset(&action.sa_mask);
 pool_signal=0;sigaction(SIGTERM,&action,&old_term);sigaction(SIGINT,&action,&old_int);
 std::vector<PoolChild> children;children.reserve(workers);const pid_t parent=getpid();
 auto restore=[&]{sigaction(SIGTERM,&old_term,nullptr);sigaction(SIGINT,&old_int,nullptr);};
 try {
  std::cout.flush();std::cerr.flush();
  for(std::size_t slot=0;slot<workers;++slot) {
   if(pool_signal)throw std::runtime_error("coordinator interrupted during startup");
   int fd[2];if(socketpair(AF_UNIX,SOCK_SEQPACKET,0,fd))throw std::runtime_error("worker socket creation failed");
   pid_t pid=fork();
   if(pid<0){close(fd[0]);close(fd[1]);throw std::runtime_error("worker fork failed");}
   if(pid==0) {
    signal(SIGTERM,SIG_DFL);signal(SIGINT,SIG_DFL);pool_signal=0;
    if(prctl(PR_SET_PDEATHSIG,SIGKILL)||getppid()!=parent)_exit(120);
    close(fd[0]);for(const auto& c:children)close(c.fd);
    cpu_set_t one;CPU_ZERO(&one);CPU_SET(cpus[slot],&one);
    if(sched_setaffinity(0,sizeof one,&one))_exit(121);
    std::ofstream log(evidence/("worker-"+std::to_string(slot)+".log"));
    log<<"pid="<<getpid()<<" cpu="<<sched_getcpu()<<" logical_affinity_count=1\n";log.flush();
    try {
     for(;;) {
      std::uint64_t index;Receive(fd[1],&index,sizeof index);
      if(index==UINT64_MAX){close(fd[1]);_exit(0);}
      try{callback(index,slot);PoolReply reply{index,1};Send(fd[1],&reply,sizeof reply);}
      catch(const std::exception& e){log<<"FAIL task="<<index<<" "<<e.what()<<'\n';log.flush();PoolReply reply{index,0};Send(fd[1],&reply,sizeof reply);_exit(1);}
     }
    }catch(const std::exception& e){log<<"FAIL IPC "<<e.what()<<'\n';log.flush();_exit(122);}
   }
   close(fd[1]);children.push_back({pid,fd[0],0,false});
  }
  std::size_t next=0,complete=0;
  auto assign=[&](std::size_t slot) {
   auto& c=children[slot];if(next>=tasks)return;
   c.task=next++;c.busy=true;Send(c.fd,&c.task,sizeof c.task);
   dispatch<<c.task<<'\t'<<slot<<'\t'<<c.pid<<'\t'<<cpus[slot]<<"\tSTART\n";dispatch.flush();
  };
  for(std::size_t i=0;i<workers;++i)assign(i);
  while(complete<tasks) {
   if(pool_signal)throw std::runtime_error("coordinator interrupted by signal "+std::to_string(pool_signal));
   std::vector<pollfd> fds;for(const auto& c:children)fds.push_back({c.fd,short(POLLIN),0});
   const int rc=poll(fds.data(),fds.size(),100);
   if(rc<0){if(errno==EINTR)continue;throw std::runtime_error("worker poll failed");}
   std::vector<std::size_t> ready;
   // Consume all ready replies before assigning more: a reported failure wins.
   for(std::size_t i=0;i<workers;++i)if(fds[i].revents) {
    auto& c=children[i];PoolReply reply{};
    Receive(c.fd,&reply,sizeof reply);
    if(!c.busy||reply.index!=c.task||reply.ok!=1)throw std::runtime_error("worker failure at task "+std::to_string(c.task));
    dispatch<<c.task<<'\t'<<i<<'\t'<<c.pid<<'\t'<<cpus[i]<<"\tCOMPLETE\n";dispatch.flush();
    c.busy=false;++complete;ready.push_back(i);
   }
   if(pool_signal)throw std::runtime_error("coordinator interrupted before dispatch");
   for(auto i:ready)assign(i);
  }
  const std::uint64_t stop=UINT64_MAX;for(auto& c:children)Send(c.fd,&stop,sizeof stop);
  for(auto& c:children){int status=0;pid_t r;do{r=waitpid(c.pid,&status,0);}while(r<0&&errno==EINTR&&!pool_signal);
   if(r!=c.pid||!WIFEXITED(status)||WEXITSTATUS(status)!=0)throw std::runtime_error("worker did not exit cleanly");
   c.pid=-1;close(c.fd);c.fd=-1;}
  restore();
 }catch(const std::exception& e) {
  std::ofstream failure(evidence/"failure.txt");failure<<e.what()<<'\n';failure.close();
  StopChildren(children);restore();throw;
 }
}
} // namespace NativeQualification

#else
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <vector>
namespace NativeQualification {
inline std::vector<int> PhysicalCpus(){throw std::runtime_error("parallel qualification requires Linux");}
inline void ProcessPool(std::size_t,std::size_t,const std::filesystem::path&,
 const std::function<void(std::size_t,std::size_t)>&){throw std::runtime_error("parallel qualification requires Linux");}
}
#endif
