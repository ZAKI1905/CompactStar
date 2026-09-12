#pragma once
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace CompactStar::Physics::Rotochemical
{
// SHA-256 is provenance only, never a numerical coefficient authority.
inline std::string SourceSHA256(const std::string& bytes)
{
    static constexpr std::array<uint32_t,64> k{
      0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
      0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
      0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
      0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
      0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
      0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
      0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
      0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    std::array<uint32_t,8> h{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    std::vector<unsigned char> b(bytes.begin(),bytes.end());const uint64_t bits=uint64_t(b.size())*8;
    b.push_back(0x80);while(b.size()%64!=56)b.push_back(0);
    for(int i=7;i>=0;--i)b.push_back(static_cast<unsigned char>(bits>>(8*i)));
    auto rot=[](uint32_t v,unsigned n){return (v>>n)|(v<<(32-n));};
    for(size_t off=0;off<b.size();off+=64){std::array<uint32_t,64>w{};
      for(size_t i=0;i<16;++i)for(size_t j=0;j<4;++j)w[i]=(w[i]<<8)|b[off+4*i+j];
      for(size_t i=16;i<64;++i){auto a=w[i-15],c=w[i-2];w[i]=w[i-16]+(rot(a,7)^rot(a,18)^(a>>3))+w[i-7]+(rot(c,17)^rot(c,19)^(c>>10));}
      auto v=h;for(size_t i=0;i<64;++i){uint32_t t1=v[7]+(rot(v[4],6)^rot(v[4],11)^rot(v[4],25))+((v[4]&v[5])^(~v[4]&v[6]))+k[i]+w[i];uint32_t t2=(rot(v[0],2)^rot(v[0],13)^rot(v[0],22))+((v[0]&v[1])^(v[0]&v[2])^(v[1]&v[2]));for(size_t j=7;j>0;--j)v[j]=v[j-1];v[4]+=t1;v[0]=t1+t2;}
      for(size_t i=0;i<8;++i)h[i]+=v[i];
    }
    std::ostringstream out;out<<std::hex<<std::setfill('0');for(auto v:h)out<<std::setw(8)<<v;return out.str();
}
struct FrozenSource
{
    const std::string path,bytes,sha256;
    explicit FrozenSource(const std::filesystem::path& p):path(std::filesystem::absolute(p).string()),bytes(Read(path)),sha256(SourceSHA256(bytes)){}
    void RequireDiskCurrent()const {if(Read(path)!=bytes)throw std::runtime_error("changed frozen source: "+path);}
    static std::string Read(const std::string& path){std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("missing source: "+path);std::string b;std::array<char,65536>buf;while(f.read(buf.data(),buf.size())||f.gcount())b.append(buf.data(),size_t(f.gcount()));if(!f.eof()||f.bad())throw std::runtime_error("source read failure: "+path);return b;}
};
}
