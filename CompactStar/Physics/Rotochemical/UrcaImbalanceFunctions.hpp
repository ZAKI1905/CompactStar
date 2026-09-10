#pragma once
#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::Rotochemical
{
// ADR-0014 §3.6; FR2005 (34)-(37), R1995 (29)-(32), u=xi/pi.
// FR2005 (37) prints pi^6 in the last HM term, inconsistent with its (59)-(60)
// and the independent Fermi convolution. Use ratified pi^8; no erratum claimed.
// Pure mathematics: no stellar, thermal, coefficient or process selection state.
struct UrcaImbalanceFunctions
{
    static double DirectIncrement(double xi)
    { const double u2 = Square(xi); return u2 * (1071 + u2 * (315 + 21*u2))/457; }
    static double ModifiedIncrement(double xi)
    { const double u2 = Square(xi); return u2 * (22020 + u2*(5670 + u2*(420 + 9*u2)))/11513; }
    static double FD(double xi) { return 1 + DirectIncrement(xi); }
    static double FM(double xi) { return 1 + ModifiedIncrement(xi); }
    static double HD(double xi)
    { const double u2=Square(xi); return xi/(Pi()*Pi())*(714+u2*(420+42*u2))/457; }
    static double HM(double xi)
    { const double u2=Square(xi); return xi/(Pi()*Pi())*(14680+u2*(7560+u2*(840+24*u2)))/11513; }
    static double Pi() { return std::acos(-1.0); }
  private:
    static double Square(double xi)
    { if(!std::isfinite(xi)) throw std::runtime_error("nonfinite Urca xi"); const double u=xi/Pi(); return u*u; }
};
}
