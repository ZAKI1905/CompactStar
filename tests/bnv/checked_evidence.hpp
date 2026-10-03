#pragma once
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <limits>
#include <cmath>

// Publish one complete TSV only after checked write/flush/close and read-back.
// The output directory itself is published by the qualification runner last.
inline void PublishTsv(const std::filesystem::path& path,const std::string& bytes,
                       std::size_t columns,std::size_t rows,int time_column=-1,
                       std::size_t fail_after=std::numeric_limits<std::size_t>::max())
{
    if(std::filesystem::exists(path))throw std::runtime_error("evidence already exists");
    const auto temp=path.string()+".partial";
    if(std::filesystem::exists(temp))throw std::runtime_error("partial evidence already exists");
    try
    {
        std::ofstream f;f.exceptions(std::ios::badbit|std::ios::failbit);f.open(temp,std::ios::binary);
        if(bytes.size()>fail_after)
        {f.write(bytes.data(),fail_after);f.setstate(std::ios::badbit);}
        f.write(bytes.data(),bytes.size());f.flush();f.close();
        std::ifstream in;in.exceptions(std::ios::badbit);in.open(temp,std::ios::binary);
        if(!in)throw std::runtime_error("evidence read-back open failure");
        std::string recovered((std::istreambuf_iterator<char>(in)),{});
        if(recovered!=bytes||bytes.empty()||bytes.back()!='\n')throw std::runtime_error("evidence byte mismatch");
        std::istringstream all(recovered);std::string line;std::size_t count=0;double last=-INFINITY;
        while(std::getline(all,line))
        {
            std::istringstream row(line);std::string cell;std::size_t n=0;
            while(std::getline(row,cell,'\t'))
            {
                if(count&&int(n)==time_column){std::size_t end;double t=std::stod(cell,&end);
                    if(end!=cell.size()||!std::isfinite(t)||!(t>last))throw std::runtime_error("evidence time order");last=t;}
                ++n;
            }
            if(n!=columns)throw std::runtime_error("evidence column count");++count;
        }
        if(count!=rows+1)throw std::runtime_error("evidence row count");
        std::filesystem::rename(temp,path);
    }
    catch(...){std::error_code ec;std::filesystem::remove(temp,ec);throw;}
}

inline void PublishTable(const std::filesystem::path& path,const std::string& bytes,
                         std::size_t rows,int time_column=-1)
{
    const auto newline=bytes.find('\n');
    if(newline==std::string::npos)throw std::runtime_error("missing table header");
    const auto columns=1+std::count(bytes.begin(),bytes.begin()+newline,'\t');
    PublishTsv(path,bytes,columns,rows,time_column);
}

inline void PublishReplacingTable(const std::filesystem::path& path,const std::string& bytes,
                                  std::size_t rows)
{
    const std::filesystem::path staged=path.string()+".replacement";
    PublishTable(staged,bytes,rows);
    // Same-directory rename replaces the old complete snapshot only after validation.
    std::filesystem::rename(staged,path);
}
