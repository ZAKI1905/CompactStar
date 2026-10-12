// Validation only: production's existing exact-profile guard, no BNV trajectory.
#include "../rotochemical/fixture.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

int main(int argc,char** argv)
{
 try {
    require(argc==4,"native-profile native-certificate fresh-work");
    const std::filesystem::path profile=argv[1],certificate=argv[2],work=argv[3];
    auto fixture=Fixture(profile,certificate.string(),work/"nominal",80000);
    fixture.g->RequireCurrent();fixture.z->RequireCurrent();
    std::cout<<"NATIVE_EXACT_PROFILE_BINDING PASS\n";
    const auto mutant=work/"mutated-profile";
    std::filesystem::copy(profile,mutant,std::filesystem::copy_options::recursive);
    std::ifstream in(mutant/"profile.tsv");std::string line,bytes;
    require(bool(std::getline(in,line)),"missing header");bytes=line+"\n";
    require(bool(std::getline(in,line)),"missing first row");
    const auto tab=line.find('\t');require(tab!=std::string::npos,"missing first column");
    const double radius=std::stod(line.substr(0,tab));
    std::ostringstream altered;altered<<std::setprecision(17)<<std::nextafter(radius,std::numeric_limits<double>::infinity())<<line.substr(tab)<<'\n';
    bytes+=altered.str();while(std::getline(in,line))bytes+=line+'\n';in.close();
    {std::ofstream out(mutant/"profile.tsv");out<<bytes;out.close();require(bool(out),"mutation write failure");}
    bool refused=false;
    try { Fixture(mutant,certificate.string(),work/"mutant",80000); }
    catch(const std::exception& e) {refused=std::string(e.what()).find("characterization profile differs from owning canonical star")!=std::string::npos;std::cout<<"EXACT_PROFILE_MUTATION_REFUSAL "<<e.what()<<'\n';}
    require(refused,"one-ULP profile mismatch survived exact guard");
    // Exercise real retained-byte currentness against a disposable certificate.
    const auto local_certificate=work/"currentness-certificate.txt";
    std::filesystem::copy_file(certificate,local_certificate);
    auto local=Fixture(profile,local_certificate.string(),work/"currentness",80000);
    local.g->RequireCurrent();local.z->RequireCurrent();
    {std::ofstream out(local_certificate,std::ios::app);out<<"\n";out.close();require(bool(out),"mutation write failure");}
    refused=false;try {local.z->RequireCurrent();}catch(const std::exception& e){refused=true;std::cout<<"CURRENTNESS_MUTATION_REFUSAL "<<e.what()<<'\n';}
    require(refused,"changed input certificate survived owner currentness");
    std::cout<<"BA13_INPUT_IDENTITY_CONTROLS PASS; no frozen sensitivity certificate attached; no trajectory\n";
    return 0;
 } catch(const std::exception& e) {std::cerr<<"STOP "<<e.what()<<'\n';return 1;}
}
