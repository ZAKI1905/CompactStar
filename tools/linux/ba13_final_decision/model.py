"""Validation-sidecar algebra. Accuracy requirements are never evidence estimates."""
import math

def correlated_increment(a, d, t, denominator):
    return (a - t*d) / denominator

def joint_ratio_bound(joint_residual, denominator_error, denominator):
    if joint_residual is None or denominator_error is None:
        raise ValueError('missing joint residual or denominator error')
    if joint_residual < 0 or denominator_error < 0 or denominator_error >= abs(denominator):
        raise ValueError('denominator not separated from zero')
    return joint_residual / (abs(denominator)-denominator_error)

def feasible_cap(threshold, drift, finite, fit_norm, base_band, *, absolute=False):
    if not all(math.isfinite(v) for v in (threshold,drift,finite,fit_norm,base_band)):
        raise ValueError('nonfinite feasibility input')
    if threshold<=0 or min(drift,finite,fit_norm,base_band)<0 or base_band>=1:
        raise ValueError('invalid feasibility input')
    slack=threshold-drift-finite
    rhs=slack if absolute else slack*(1-base_band)
    coefficient=2+fit_norm+(0 if absolute else threshold)
    return dict(pointwise_coefficient=coefficient,rhs=rhs,
                single_axis_pointwise_cap=rhs/coefficient,
                single_axis_gap_or_fit_cap=rhs,
                evidence_type='required_accuracy_not_estimate',
                joint_constraint='coefficient*e + g + z <= rhs',
                qualifying_value=None)

def continuous_remainder(curvature, width, *, reference_kind, evidence_scope):
    if reference_kind != 'X2':
        raise ValueError('wrong reference attachment')
    if evidence_scope not in ('analytic_supremum','accepted_conditional_supremum'):
        raise ValueError('sampled curvature or holdouts alone are insufficient')
    if curvature is None or width < 0 or curvature < 0:
        raise ValueError('missing or invalid curvature')
    return curvature*width*width/8

def complete_total(terms, registry):
    values=[]
    for t in terms:
        v=t.get('value'); kind=t.get('evidence_type')
        if v is None or kind=='unresolved' or t.get('missing'):
            return None
        if not math.isfinite(v) or v<0 or not t.get('provenance'):
            raise ValueError('invalid contribution')
        if kind in ('conditional_estimate','validated_numerical_estimate'):
            field='condition_id' if kind=='conditional_estimate' else 'standard_id'
            table='accepted_conditions' if field=='condition_id' else 'accepted_standards'
            authority=t.get('authority')
            if not authority or registry.get(table,{}).get(t.get(field))!=authority:
                return None
        elif kind=='analytic_bound':
            if not t.get('derivation'):raise ValueError('missing derivation')
        elif kind!='measured_difference':
            raise ValueError('invalid evidence class')
        values.append(v)
    names={t['name'] for t in terms}
    if any(set(t.get('dependencies',()))-names for t in terms):return None
    return math.fsum(values)
