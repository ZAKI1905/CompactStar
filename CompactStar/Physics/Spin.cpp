// These declared diagnostics had no implementation in the historical archive.
// Python's dynamic_lookup concealed the missing symbols at executable link time.
// Keep them explicitly unavailable until their scientific contract is qualified.
#include "CompactStar/Physics/Spin.hpp"
#include <stdexcept>

namespace CompactStar::Physics::Spin
{
double CharacteristicAge(const State::SpinState &)
{
    throw std::logic_error("CharacteristicAge: diagnostic implementation not qualified");
}

double DipoleFieldEstimate(const State::SpinState &, Core::StarProfileView)
{
    throw std::logic_error("DipoleFieldEstimate: normalization not specified or qualified");
}
}
