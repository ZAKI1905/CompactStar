#pragma once
#include <array>
#include <stdexcept>
#include <string>

namespace CompactStar::Physics::Rotochemical
{
// Byte identities, not scientific certificates. Selection never changes the
// owning-star, lifetime, revision, disk-currentness or BNV certificate checks.
// No caller-supplied replacement hashes are accepted.
enum class QualifiedBackgroundInputs { HistoricalStructure1, NativeStructure1 };
inline std::array<std::string,4> BackgroundInputHashes(QualifiedBackgroundInputs identity)
{
    switch(identity) {
    case QualifiedBackgroundInputs::HistoricalStructure1:
        return {{"e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890",
                 "3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4",
                 "7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453",
                 "7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a"}};
    case QualifiedBackgroundInputs::NativeStructure1:
        return {{"62b8dc594fd9c30d413ccb1d906459e2599049966b0bec097284b43598a07a56",
                 "6dc75c841408d8e03fed5cbc8444310af35c77ab34d0339df309c00cbb3e7f4a",
                 "337ed8ae8d6168d8a612f1344f204dc927ac9e0867a47542d5b943dbf27b8537",
                 "f6b9c4068a430cc356b54e39c21a0ca17e17b8c2b35d17b836638298d1112d7f"}};
    }
    throw std::runtime_error("unknown qualified background input identity");
}
}
