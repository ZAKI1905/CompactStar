#include "checkpoint_transport.hpp"
#include <cstring>
#include <iostream>
#include <sstream>
#include <limits>
using namespace NativeQualification::Transport;
struct Fill {
 std::size_t counter=0;
 template<class T> void Set(T& v) {
  if constexpr(std::is_same_v<T,double>) {
   const double values[]={-0.0,std::numeric_limits<double>::denorm_min(),std::nextafter(1.,2.),1e-27,-1e42};
   v=values[counter++%5];
  }else if constexpr(std::is_arithmetic_v<T>||std::is_enum_v<T>)v=static_cast<T>(1);
  else Fields(*this,v);
 }
 void Set(std::string& v){v="exact IPC field "+std::to_string(++counter);}
 template<class T,std::size_t N> void Set(std::array<T,N>& v){for(auto& x:v)Set(x);}
 template<class T> void Set(std::vector<T>& v){v.resize(2);for(auto& x:v)Set(x);}
 template<class K,class V> void Set(std::map<K,V>& v){K k;V x;Set(k);Set(x);v.emplace(k,x);}
 template<class... T> void operator()(T&... v){(Set(v),...);}
};
void Check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(int argc,char** argv) {
 try {
  Check(argc==2,"fresh output required");const std::filesystem::path root=argv[1];
  Check(!std::filesystem::exists(root),"fresh output required");std::filesystem::create_directories(root);
  B::CheckpointOutput q;Fill fill;fill(q);Write(root/"a.bin",q);const auto r=Read(root/"a.bin");Write(root/"b.bin",r);
  auto bytes=[](const auto& p){std::ifstream f(p,std::ios::binary);return std::string(std::istreambuf_iterator<char>(f),{});};
  Check(bytes(root/"a.bin")==bytes(root/"b.bin"),"binary roundtrip mismatch");
  Check(std::memcmp(q.state.data(),r.state.data(),sizeof(q.state))==0,"binary64 bits changed");
  for(const auto& mutation:{std::string("truncated"),std::string("trailing"),std::string("schema")}) {
   auto data=bytes(root/"a.bin");if(mutation=="truncated")data.resize(data.size()/2);
   if(mutation=="trailing")data+='!';if(mutation=="schema")data.at(8)='!';
   {std::ofstream f(root/mutation,std::ios::binary);f.write(data.data(),data.size());}
   bool refused=false;try{Read(root/mutation);}catch(const std::exception&){refused=true;}Check(refused,"malformed packet accepted");
  }
  std::cout<<"PASS all IPC fields, exact binary64 including signed zero/subnormal, truncation/trailing/schema refusals\n";return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
