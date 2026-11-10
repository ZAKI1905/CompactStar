#pragma once
#include <cmath>
#include <stdexcept>

namespace NativeQualification {
// Test-harness policy only. The immutable accepted history is never altered.
// Choose the latest positive-prefix endpoint leaving at least half the horizon.
template<class History>
typename History::const_iterator RestartAtOrBeforeMidpoint(
    const History& history, double start, double end) {
    if(!std::isfinite(start)||!std::isfinite(end)||!(end>start))
        throw std::runtime_error("invalid restart horizon");
    const double midpoint=start+(end-start)/2;
    if(!std::isfinite(midpoint)||!(midpoint>start&&midpoint<end))
        throw std::runtime_error("unrepresentable restart midpoint");
    auto selected=history.cend();
    double previous=start;
    for(auto it=history.cbegin();it!=history.cend();++it) {
        if(it->t_left_s!=previous||!std::isfinite(it->t_right_s)
           ||!(it->t_right_s>previous)||it->t_right_s>end)
            throw std::runtime_error("invalid accepted restart history");
        if(it->t_right_s<=midpoint)selected=it;
        previous=it->t_right_s;
    }
    if(previous!=end||selected==history.cend())
        throw std::runtime_error("no accepted restart endpoint at or before midpoint");
    return selected;
}
} // namespace NativeQualification
