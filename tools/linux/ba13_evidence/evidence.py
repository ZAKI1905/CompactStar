"""Strict successor decision layer; historical classifiers are replay-only."""
import math

FRACTIONS = dict(eos_equilibrium=.08,center=.04,surface_events=.06,integration=.10,
 implicit_stages=.02,floating_point=.02,common_B=.10,production_representation=.23,
 differentiation_inverse=.10,normalization=.05,nonlinear_fit=.10,between_samples=.10)
CLASSES = {'analytic_bound','measured_difference','validated_numerical_estimate',
           'conditional_estimate','unresolved'}


def contribution(name, category, value, *, evidence_type, provenance,
                 condition_id=None, standard_id=None, authority=None,
                 derivation=None, dependencies=(), missing=None):
    if evidence_type not in CLASSES or category not in FRACTIONS:
        raise ValueError('explicit recognized class/category required')
    if value is not None and (not math.isfinite(value) or value < 0):
        raise ValueError('nonfinite/negative uncertainty')
    if not provenance or not all(p.get('path') and len(p.get('sha256',''))==64 for p in provenance):
        raise ValueError('source path and SHA256 provenance required')
    if evidence_type == 'analytic_bound' and not derivation:
        raise ValueError('analytic evidence needs a derivation and explicit input dependencies')
    if evidence_type == 'conditional_estimate' and not condition_id:
        raise ValueError('conditional evidence needs an itemized condition')
    if evidence_type == 'validated_numerical_estimate' and not standard_id:
        raise ValueError('numerical evidence needs a named standard')
    if value is None and not missing:
        raise ValueError('missing term needs a reason')
    return dict(name=name,category=category,value=value,evidence_type=evidence_type,
                provenance=provenance,condition_id=condition_id,standard_id=standard_id,
                authority=authority,derivation=derivation,dependencies=list(dependencies),missing=missing)


def qualify(history, terms, registry, *, frozen_historical, exact=False):
    if history != frozen_historical:
        raise ValueError('historical identity changed')
    if set(t['category'] for t in terms) != set(FRACTIONS):
        raise ValueError('missing contribution category')
    if len({t['name'] for t in terms}) != len(terms):
        raise ValueError('duplicate contribution name')
    by_name={t['name']:t for t in terms}
    reasons={}
    for t in terms:
        # Revalidate serialized records as well as records from the constructor.
        contribution(**t)
        n=t['name']; kind=t['evidence_type']
        if t['value'] is None or kind=='unresolved' or t['missing']:
            reasons[n]=t['missing'] or 'unresolved evidence'
        elif kind in ('conditional_estimate','validated_numerical_estimate'):
            field='condition_id' if kind=='conditional_estimate' else 'standard_id'
            table='accepted_conditions' if field=='condition_id' else 'accepted_standards'
            acceptance=registry.get(table,{}).get(t[field])
            if not t['authority'] or not acceptance or acceptance!=t['authority']:
                reasons[n]='missing recorded acceptance: '+str(t[field])
        for dep in t['dependencies']:
            if dep not in by_name:
                reasons[n]='missing dependency: '+dep
    # Propagation cannot upgrade its inputs, including transitive dependencies.
    for _ in terms:
        for t in terms:
            if any(d in reasons for d in t['dependencies']):
                reasons[t['name']]='unqualified input dependency'
    added=sum(t['value'] or 0. for t in terms)
    h=history['threshold']-history['drift_bound']
    category={c:sum(t['value'] or 0. for t in terms if t['category']==c) for c in FRACTIONS}
    measured={c:sum(t['value'] or 0. for t in terms if t['category']==c and
                   t['evidence_type']=='measured_difference') for c in FRACTIONS}
    violated=[c for c,v in category.items() if v>h*FRACTIONS[c]]
    measured_failure=any(v>h*FRACTIONS[c] for c,v in measured.items())
    if exact and (history['drift_bound']!=0 or any(t['value'] not in (None,0) for t in terms)):
        measured_failure=True
    if exact and any(t['evidence_type']!='analytic_bound' for t in terms):
        reasons['exact']='exact identity cannot use estimates'
    status='FAIL' if measured_failure else 'UNQUALIFIED' if reasons else (
        'FAIL' if violated or added>h else 'PASS')
    return dict(status=status,refusal_reasons=reasons,conditional_terms=[t['name'] for t in terms
                if t['evidence_type']=='conditional_estimate'],missing_terms=[t['name'] for t in terms if t['value'] is None],
                finite_added_estimates=added,added_justified_error=None if reasons else added,
                partial_utilization=(history['drift_bound']+added)/history['threshold'],
                proposed_total_utilization=None if reasons else (history['drift_bound']+added)/history['threshold'],
                remaining_physical_headroom=None if reasons else h-added,
                remaining_headroom_after_finite_terms=h-added,category_totals=category,
                allocated_budget={c:h*f for c,f in FRACTIONS.items()},category_budget_exceeded=violated)


def m1_decision(discrepancy, reference_band, share, *, complete, accepted):
    if not math.isfinite(discrepancy) or discrepancy<0 or not math.isfinite(share) or share<0:
        raise ValueError('invalid comparison')
    if discrepancy>share:
        return 'FAIL'
    if reference_band is None or not complete or not accepted:
        return 'UNQUALIFIED'
    if not math.isfinite(reference_band) or reference_band<0:
        raise ValueError('invalid reference band')
    return 'PASS' if discrepancy+reference_band<=share else 'UNQUALIFIED'
