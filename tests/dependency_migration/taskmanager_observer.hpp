// Test-only read access; never instantiated in the production library.
#include <Confind/ContourFinder.hpp>
#include <Zaki/Math/Math_Core.hpp>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <stdexcept>
namespace CSMObserve {
inline FILE* output() {
    static FILE* f = [] {
        const char* name = std::getenv("CSM_OBSERVATIONS");
        FILE* result = name ? std::fopen(name, "w") : nullptr;
        if (!result) throw std::runtime_error("CSM_OBSERVATIONS required");
        return result;
    }();
    return f;
}
inline void value(const char* key, double v) {
    std::uint64_t bits; std::memcpy(&bits, &v, sizeof(bits));
    std::fprintf(output(), "%s\t%a\t%016llx\n", key, v, static_cast<unsigned long long>(bits));
    std::fflush(output());
}
inline void point(const char* key, const Zaki::Math::Coord2D& p) {
    value(key, p.x); value(key, p.y);
}
inline void curve(const char* key, const Zaki::Math::Curve2D& c) {
    value(key, c.Size()); for (const auto& p : c.pts) point(key, p);
}
inline void points(const char* key, const std::vector<Zaki::Math::Coord2D>& ps) {
    value(key, ps.size()); for (const auto& p : ps) point(key, p);
}
inline void contours(const char* key, const std::vector<CONFIND::Cont2D>& cs) {
    value(key, cs.size());
    for (const auto& c : cs) {
        value(key, c.GetVal()); value(key, c.GetFound()); value(key, c.size());
        for (size_t i=0; i<c.size(); ++i) {
            const auto p = c[i]; value(key, p.x); value(key, p.y); value(key, p.z);
        }
    }
}
}
