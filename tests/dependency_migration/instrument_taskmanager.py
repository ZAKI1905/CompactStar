#!/usr/bin/env python3
"""Generate a test-only observer TU from actual TaskManager source.

No arithmetic expression or control condition is replaced. Fail on anchor drift.
The ordinary integration executable remains linked to the uninstrumented library.
"""
import sys
from pathlib import Path
source, output = map(Path, sys.argv[1:])
s = source.read_text()
def after(anchor, observation, count=1):
    global s
    if s.count(anchor) != count:
        raise RuntimeError(f'observation anchor drift: {anchor}')
    s = s.replace(anchor, anchor + '\n' + observation)
after('std::vector<CONFIND::Cont2D> con_set = con.GetContourSet();', 'CSMObserve::contours("critical.raw", con_set);')
after('critical_curve[-1] = seg.P(lam);', 'CSMObserve::value("critical.argmax.contour", c); CSMObserve::value("critical.argmax.point", p); CSMObserve::value("critical.argmax.mass", tmp_max_mass_tot); CSMObserve::point("critical.selected", critical_curve[-1]);')
after('critical_curve.MakeSmooth(10);', 'CSMObserve::curve("critical.smooth", critical_curve);')
after('critical_curve.Append({sequence_grid[2][m_max_v_idx], sequence_grid[8][m_max_v_idx]});', 'CSMObserve::curve("critical.final", critical_curve);')
after('std::vector<CONFIND::Cont2D> mass_con_set = con.GetContourSet();', 'CSMObserve::contours("mass.raw", mass_con_set);')
after('Zaki::Math::Curve2D mass_curve = mass_con_set[0].ConvertToCurve2D();', 'CSMObserve::curve("mass.curve", mass_curve); CSMObserve::point("mass.first", mass_curve[0]);')
after('std::vector<Zaki::Math::Coord2D> intersection = critical_curve.Intersection(mass_curve);', 'CSMObserve::curve("critical.reimported", critical_curve); CSMObserve::points("mass.intersection", intersection);')
after('double cont_lvl_max = B_tot_grid.Evaluate(in_m_range[1].x, in_m_range[1].y);', 'CSMObserve::value("baryon.min", cont_lvl_min); CSMObserve::value("baryon.max", cont_lvl_max);')
after('std::vector<CONFIND::Cont2D> B_con_set = con.GetContourSet();', 'CSMObserve::contours("baryon.raw", B_con_set);')
after('Zaki::Math::Curve2D B_curve = B_con_set[i].ConvertToCurve2D();', 'CSMObserve::curve("baryon.curve", B_curve);')
after('std::vector<Zaki::Math::Coord2D> intersection = critical_curve.Intersection(B_curve);', 'CSMObserve::points("baryon.intersection", intersection);')
after('B_curv_stable.emplace_back(B_curve.Bisect(intersection[0]).first);', 'CSMObserve::value("baryon.getidx", B_curve.GetIdx(intersection[0])); CSMObserve::curve("baryon.bisect", B_curv_stable.back());')
s = '#include "taskmanager_observer.hpp"\n' + s
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(s)
