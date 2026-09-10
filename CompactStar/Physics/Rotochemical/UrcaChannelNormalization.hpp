#pragma once
#include <CompactStar/Physics/Rotochemical/ChemicalImbalanceState.hpp>
#include <array>
#include <functional>
#include <initializer_list>
#include <string>
#include <vector>

namespace CompactStar::Physics::Rotochemical
{
enum class UrcaProcess { De, Dmu, Me, Mmu };
inline std::size_t ProcessIndex(UrcaProcess p)
{ switch(p) {case UrcaProcess::De:return 0;case UrcaProcess::Dmu:return 1;case UrcaProcess::Me:return 2;case UrcaProcess::Mmu:return 3;} throw std::runtime_error("invalid Urca process"); }
inline unsigned TemperatureExponent(UrcaProcess p) { return ProcessIndex(p)<2?6:8; }
inline BetaChannel LeptonChannel(UrcaProcess p) { return ProcessIndex(p)%2==0?BetaChannel::Npe:BetaChannel::NpMu; }
inline constexpr std::array<UrcaProcess,4> UrcaProcesses{UrcaProcess::De,UrcaProcess::Dmu,UrcaProcess::Me,UrcaProcess::Mmu};
class UrcaProcessSelection
{
  public:
    explicit UrcaProcessSelection(std::initializer_list<UrcaProcess> enabled)
    { for(auto p:enabled) { auto i=ProcessIndex(p); if(enabled_[i]) throw std::runtime_error("duplicate Urca process"); enabled_[i]=true; } }
    static UrcaProcessSelection ControlledModifiedOnly() { return UrcaProcessSelection{UrcaProcess::Me,UrcaProcess::Mmu}; }
    bool Enabled(UrcaProcess p) const { return enabled_[ProcessIndex(p)]; }
  private:
    std::array<bool,4> enabled_{};
};
struct UrcaSupportInterval { double left_km, right_km; };
// Applicability is supplied separately. A true triangle never enables a process.
inline bool DirectUrcaTriangle(double nn,double np,double nl)
{
    if(!(nn>0 && np>0 && nl>0)) return false;
    return std::cbrt(nn)<=std::cbrt(np)+std::cbrt(nl);
}
// Layer 1: local normalization only. Provider functions are sampled at construction
// of the global result; no evolving support/temperature predicate is accepted.
class UrcaChannelNormalization
{
  public:
    using LocalFunction = std::function<double(double)>;
    UrcaChannelNormalization(UrcaProcess process, LocalFunction source,
        std::vector<UrcaSupportInterval> support, std::string identity)
      : process_(process),source_(std::move(source)),support_(std::move(support)),identity_(std::move(identity))
    {
        ProcessIndex(process_);
        if(!source_ || identity_.empty()) throw std::runtime_error("normalization provenance required");
        double last=-1;
        for(auto s:support_) {
            if(!std::isfinite(s.left_km)||!std::isfinite(s.right_km)||s.left_km<0||s.right_km<=s.left_km||s.left_km<last)
                throw std::runtime_error("support must be disjoint innermost-first intervals");
            last=s.right_km;
        }
    }
    UrcaProcess Process() const { return process_; }
    const std::vector<UrcaSupportInterval>& Support() const { return support_; }
    const std::string& Identity() const { return identity_; }
    // Units erg cm^-3 s^-1 K^-q. Positive on declared enabled support.
    double LocalS(double r_km) const
    { double s=source_(r_km); if(!(s>0)||!std::isfinite(s)) throw std::runtime_error("nonpositive local Urca normalization"); return s; }
  private:
    UrcaProcess process_;
    LocalFunction source_;
    std::vector<UrcaSupportInterval> support_;
    std::string identity_;
};
}
