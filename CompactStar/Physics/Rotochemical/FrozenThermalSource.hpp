#pragma once
#include <CompactStar/Physics/Rotochemical/FrozenSource.hpp>
#include <CompactStar/EOS/CompOSE_Thermo.hpp>
#include <memory>
#include <functional>
namespace CompactStar::Physics::Rotochemical
{
// Const allocation prevents a caller from retaining a mutable parser alias.
class FrozenThermalSource final
{
    friend struct FrozenThermalTestAccess;
  public:
    FrozenThermalSource(const std::filesystem::path& directory,EOS::CompOSE_Thermo::Options options,std::string identity):FrozenThermalSource(directory,options,std::move(identity),{}){}
  private:
    FrozenThermalSource(const std::filesystem::path& directory,EOS::CompOSE_Thermo::Options options,std::string identity,const std::function<void()>& after_parse):identity_(std::move(identity))
    {
        if(identity_.empty()||options.clamp_to_domain||options.Tmin_for_derivative_MeV!=0)throw std::runtime_error("controlled thermal identity/options mismatch");
        for(const char* n:{"eos.t","eos.nb","eos.yq","eos.thermo"})sources_.emplace_back(directory/n);
        table_=std::make_shared<const EOS::CompOSE_Thermo>(directory.string(),options);
        if(!table_->IsLoaded())throw std::runtime_error("thermal parse failed");
        if(table_->TGrid_MeV()!=std::vector<double>{0,1,2,4}||table_->YqGrid()!=std::vector<double>{0,1})throw std::runtime_error("controlled thermal axes mismatch");
        for(double n:table_->NbGrid_fm3())if(!(n>0)||!std::isfinite(n))throw std::runtime_error("invalid thermal density axis");
        if(after_parse)after_parse();
        RequireDiskCurrent();
    }
  public:
    void RequireDiskCurrent()const {for(const auto& s:sources_)s.RequireDiskCurrent();}
    const EOS::CompOSE_Thermo& Table()const noexcept{return *table_;}
    const std::string& Identity()const noexcept{return identity_;}
    const std::vector<FrozenSource>& Sources()const noexcept{return sources_;}
  private:
    const std::string identity_;
    std::vector<FrozenSource> sources_;
    std::shared_ptr<const EOS::CompOSE_Thermo> table_;
};
}
