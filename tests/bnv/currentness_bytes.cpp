#include <CompactStar/Analysis/ParticleNumberResponse.hpp>
#include <CompactStar/Core/NStar.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
int main(){try{
    const auto dir=std::filesystem::temp_directory_path()/"compactstar-currentness-block-read";
    if(!std::filesystem::create_directory(dir))throw std::runtime_error("test directory collision");
    struct Cleanup{std::filesystem::path p;~Cleanup(){std::error_code e;std::filesystem::remove_all(p,e);}} cleanup{dir};
    CompactStar::Core::NStar star;
    auto eos=std::make_shared<CompactStar::Analysis::NumberEosSource>();
    eos->identity="binary-currentness-test";eos->revision="1";eos->physical_domain="test";eos->table_path=(dir/"source.bin").string();
    CompactStar::Analysis::NumberProvenance p;
    p.star=&star;p.profile=&star.Profile();p.profile_version=star.Profile().Version();p.eos=eos;p.eos_snapshot=*eos;
    for(std::size_t size:{0u,1u,65535u,65536u,65537u,200003u}){
      std::string bytes(size,'\0');for(std::size_t i=0;i<size;++i)bytes[i]=static_cast<char>(i%256);
      {std::ofstream out(eos->table_path,std::ios::binary);out.exceptions(std::ios::failbit|std::ios::badbit);out.write(bytes.data(),bytes.size());out.close();}
      p.table_contents=bytes;p.RequireCurrent();
      if(size){auto changed=bytes;changed[size/2]^=1;
        {std::ofstream out(eos->table_path,std::ios::binary);out.write(changed.data(),changed.size());}
        bool rejected=false;try{p.RequireCurrent();}catch(...){rejected=true;}
        if(!rejected)throw std::runtime_error("changed source bytes escaped");}
    }
    std::filesystem::remove(eos->table_path);bool missing=false;try{p.RequireCurrent();}catch(...){missing=true;}
    if(!missing)throw std::runtime_error("missing source escaped");
    std::cout<<"PASS full-byte currentness: empty, binary, exact/adjacent block sizes, mutation, missing file\n";
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
