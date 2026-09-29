#include "CompactStar/Physics/Spin.hpp"
#include <stdexcept>
int main() {
    CompactStar::Physics::State::SpinState state;
    int rejected = 0;
    try { (void)CompactStar::Physics::Spin::CharacteristicAge(state); }
    catch (const std::logic_error&) { ++rejected; }
    try { (void)CompactStar::Physics::Spin::DipoleFieldEstimate(state, {}); }
    catch (const std::logic_error&) { ++rejected; }
    return rejected == 2 ? 0 : 1;
}
