// Replay exact binary64 inputs captured from the real historical TaskManager.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <Zaki/Math/Math_Core.hpp>
using namespace Zaki::Math;
static double read_double(std::istream& in) {
    std::uint64_t bits;in>>std::hex>>bits>>std::dec;
    if(!in)throw std::runtime_error("invalid fixture");double d;std::memcpy(&d,&bits,8);return d;
}
static Curve2D read_curve(std::istream& in) {
    size_t count;in>>count;if(!in||count==0||count>100000)throw std::runtime_error("invalid curve");
    Curve2D c;c.Reserve(count);
    for(size_t i=0;i<count;++i){double x=read_double(in),y=read_double(in);c.Append({x,y});}
    return c;
}
static FILE* out;
static size_t ordinal;
static void emit(const char* section,double d) {
    std::uint64_t bits;std::memcpy(&bits,&d,8);
    std::fprintf(out,"%s\t%zu\t%a\t%016llx\n",section,ordinal++,d,(unsigned long long)bits);
}
int main(int argc,char** argv) {
    if(argc!=3)return 2;std::ifstream input(argv[1]);if(!input)return 3;
    out=std::fopen(argv[2],"w");if(!out)return 4;
    char kind;size_t intersections=0,queries=0,levels=0;
    while(input>>kind) {
        if(kind=='L') {
            size_t n;input>>n;for(size_t i=0;i<n;++i)emit("level",read_double(input));++levels;
        } else if(kind=='I') {
            auto a=read_curve(input),b=read_curve(input);auto points=a.Intersection(b);
            emit("intersection_count",points.size());
            for(const auto& p:points){emit("intersection_x",p.x);emit("intersection_y",p.y);}++intersections;
        } else if(kind=='G') {
            auto curve=read_curve(input);double x=read_double(input),y=read_double(input);Coord2D p{x,y};
            emit("query_x",x);emit("query_y",y);emit("getidx",curve.GetIdx(p));
            auto parts=curve.Bisect(p);emit("bisect_first_size",parts.first.Size());emit("bisect_second_size",parts.second.Size());
            for(const auto& q:parts.first.pts){emit("first_x",q.x);emit("first_y",q.y);}
            for(const auto& q:parts.second.pts){emit("second_x",q.x);emit("second_y",q.y);}++queries;
        } else throw std::runtime_error("unknown fixture record");
    }
    std::fclose(out);
    return intersections==11&&queries==10&&levels==3?0:5;
}
