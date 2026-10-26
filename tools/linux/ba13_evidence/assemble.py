#!/usr/bin/env python3
"""Static successor sidecar. Finite observations never fill unknown remainders."""
import csv
import json
import os
import pathlib
import sys
import numpy as np
from reference import (S,R,O,F,RAW,LD,read,sha,arr,load_reference,interpolation,paired_comparison)
from evidence import FRACTIONS,contribution,qualify,m1_decision

sys.path.insert(0,str(S/'tools/linux'))
from ba13_refit import load


def projection_all(xs,values,errors):
    """21 evaluation rows, 11 fitted rows, unchanged historical weights."""
    xs=np.asarray(xs,float)
    sig=np.maximum(errors[10:],64*np.finfo(float).eps*np.maximum(1.,abs(values[10:])))
    w=(sig.min()/sig)**2
    z=(xs-xs[10:].mean())/np.ptp(xs[10:])
    A=np.column_stack([np.ones(21),z]);Ao=A[10:]
    H=A@np.linalg.solve(Ao.T@(w[:,None]*Ao),Ao.T*w)
    M=-H.copy();M[np.arange(10,21),np.arange(11)]+=1
    norm=max(np.max(abs(M[10:]).sum(axis=1)),np.max(1+abs(H[:10]).sum(axis=1)))
    return H,float(norm)


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Slurm required')
    comp,m1,m3,m4,diag,x,refs,K=load_reference()
    cert=read(R/'ba13-recovery/run1/fit/certificate.json')
    historical=read(F/'matrix.json')
    oldrows={(r['index'],r['quantity']):r for r in historical['rows']}
    registry=read(S/'docs/validation/ba13-evidence/acceptance-registry.json')
    tangents=read(O/'tangents.json');holdouts=read(O/'holdouts.json')
    if holdouts['status']!='PASS_EXPERIMENT':
        raise RuntimeError('reference transfer stopped; partial holdout report must be used')
    _,series=load(RAW)
    B0=LD(m4['production_B'][0]);knots=x[[0,10,20]]
    # Scalar identity evaluated on the actual nodewise historical scales.
    mu_scale=LD(comp['original_units']['mu_B']['scales'][0])
    profile_scale=arr(comp['original_units']['mu_n_profile']['scales'])
    refs['mu_n_profile']=refs['mu_B']*mu_scale/profile_scale[None,:]
    K['mu_n_profile']=K['mu_B']*float(np.max(mu_scale/profile_scale))
    rows=[];details={};m1_new={}
    source_files=[F/'component-data.json.gz',F/'m1.json',F/'m3.json',F/'m4.json',
        O/'tangents.json',O/'holdouts.json',O/'frozen-holdouts.json',
        S/'docs/validation/ba13-evidence/PREDECLARATION.md',
        S/'docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md',
        S/'tests/bnv/compare_candidate.py',
        S/'docs/validation/x2-fixed-step/evidence/ablations.json']
    provenance=[dict(path=str(p.relative_to(R)),sha256=sha(p)) for p in source_files]
    def term(name,category,value,kind,**kwargs):
        return contribution(name,category,value,evidence_type=kind,provenance=provenance,**kwargs)
    def zero(name,category,why):
        return term(name,category,0.,'analytic_bound',derivation=why)
    def missing(name,category,why):
        return term(name,category,None,'unresolved',missing=why)
    def conditional(name,category,value,condition,dependencies=()):
        return term(name,category,value,'conditional_estimate',condition_id=condition,dependencies=dependencies)
    def numerical(name,category,value,standard):
        return term(name,category,value,'validated_numerical_estimate',standard_id=standard)
    for n in series:
        vals,errs,scales,kind=series[n]
        P=np.asarray(comp['P'][n],dtype=LD)
        vals=np.asarray(vals,float);errs=np.asarray(errs,float)
        if n in refs:
            base='mu_B' if n=='mu_n_profile' else n
            mu_factor=float(np.max(mu_scale/profile_scale)) if n=='mu_n_profile' else 1.
            IX,_,_=interpolation(x,knots,refs[n],reference_kind='X2')
            slopes=np.asarray(m4['common_B'][base]['sensitivities'],dtype=LD)*mu_factor
            maxslope=float(np.max(abs(slopes)))
            root_each={str(t):maxslope*(m1['root_alignment'][str(t)]['computed_B_residual']+
                m1['root_alignment'][str(t)]['functional_GL4_tail_estimate'])/float(B0)/.999
                       for t in (0,10,20)}
            root_pair=2*max(root_each.values())
            E,D,root_pair=paired_comparison(P,IX,root_pair,reference_kind='X2')
            curve=m3['curvatures'][base]
            all_curv=[float(np.max(abs(np.asarray(w['curvature']))))*mu_factor
                      for d in curve['width_checks'].values() for w in d]
            all_curv.append(diag['B_curvature_and_interpolated_reference'][base]['target_curvature_max']*mu_factor)
            continuation_checks=bool(min(all_curv)>0 and max(all_curv)/min(all_curv)<=2 and
                all(min(m1['root_alignment'][str(t)]['slope_secants'])>0 for t in (0,10,20)))
            interps=[]
            for j in range(21):
                a,b=(0,10) if j<=10 else (10,20)
                interps.append(K[n]*float(abs((x[j]-x[a])*(x[j]-x[b])))/2)
            int_max=max(interps)
            if n.startswith('t_'):
                tail=max(float(tangents['records'][f'{t}:{n}']['fine_tail']) for t in (0,10,20))
            else:
                tail=max(m1['reference'][str(t)][base]['integration'] for t in (0,10,20))*mu_factor
            eos=max(m1['reference'][str(t)][base]['eos_equilibrium'] for t in (0,10,20))*mu_factor
            center=max(m1['reference'][str(t)][base]['center'] for t in (0,10,20))*mu_factor
            stencil=max(m1['reference'][str(t)][base]['stencil'] for t in (0,10,20))*mu_factor
            det_fit=norm=0.
            for c in range(P.shape[1]):
                H,vnorm=projection_all(x,vals[:,c],errs[:,c])
                err=np.asarray(E[:,c],float)
                det_fit=max(det_fit,float(np.max(abs(err-H@err[10:]))))
                norm=max(norm,vnorm)
            known_band=2*(tail+eos+center+stencil)+root_pair+int_max
            fit_uncertain=norm*known_band
            oldfit=m4['fit'][base]
            fit_design=(oldfit['design_perturbation']+oldfit['weight_floor_perturbation']+
                        oldfit['inner_half_excess'])*mu_factor
            u0=float(np.max(abs(P[0]-refs[n][0])))+tail+eos+center+stencil+max(root_each.values())
            if u0>=1:
                raise RuntimeError('normalization denominator not separated from zero')
            details[n]=dict(reference='X2',reference_values=refs[n].tolist(),interpolated=IX.tolist(),
                signed_error=E.tolist(),paired=D.tolist(),K=K[n],interpolation=interps,
                root_each=root_each,root_pair=root_pair,fit_operator_norm=norm,
                fit_deterministic=det_fit,fit_input_estimate=fit_uncertain,fit_design_observation=fit_design,
                base_error_estimate=u0,curvature_ratio=max(all_curv)/min(all_curv),
                continuation_sample_checks=continuation_checks,validated_supremum=False,
                reference_unseen_bias=None,
                common_B_attribution='P-I_X is at achieved B_P; reconstruction displacement c is diagnostic only; X2 residual plus B-functional tail is charged',
                source_fingerprint=sha(F/'component-data.json.gz'))
            if not continuation_checks:
                details[n]['stop']='continuation sample condition failed; no additional numerical work'
        for i in range(21):
            h={k:cert['samples'][i][k][n] for k in ('drift_bound','threshold','utilization')}
            terms=[]
            if n in refs:
                terms.append(term('direct_P_minus_interpolated_X2','production_representation',float(D[i]),'measured_difference'))
                terms.append(conditional('reference_interpolation','production_representation',max(interps[:i+1]),
                                         'global-curvature-and-cut-law'))
                terms.append(numerical('X2_discretization','integration',2*tail if i else 0.,
                    'offset-H-GL4-three-resolution-v1' if n.startswith('t_') else 'retained-GL4-three-resolution-v1'))
                terms.append(numerical('X2_EOS_sensitivity','eos_equilibrium',2*eos if i else 0.,'retained-finite-sensitivities-v1'))
                terms.append(numerical('X2_center_sensitivity','center',2*center if i else 0.,'retained-finite-sensitivities-v1'))
                terms.append(conditional('X2_common_B_alignment','common_B',root_pair if i else 0.,
                                         'positive-B-slope-and-reference-B-remainder'))
                terms.append(conditional('tangent_stencil','differentiation_inverse',2*stencil if i else 0.,
                                         'uniform-tangent-width-remainder') if n.startswith('t_') else
                             zero('direct_nonlinear_quantity','differentiation_inverse','X2 evaluates full Z inverse before comparison; no M5 duplicate'))
                terms.append(missing('reference_stage_remainder','implicit_stages',
                    'retained binary64/Newton ablations are finite primary-output sensitivities, not quantity-specific uniform stage enclosures'))
                terms.append(missing('reference_arithmetic_remainder','floating_point',
                    'primary-output precision/summation ablations do not bound every passive functional or differentiated cancellation'))
                terms.append(zero('production_events_in_comparison','surface_events',
                    'all production event/selection errors at 21 samples are in P-I_X; X2 terminates at prescribed Hs'))
                terms.append(missing('moving_cut_input_remainder','surface_events',
                    'retained cut-law residual and pressure conversion are observations; continuous law/response needs itemized acceptance'))
                terms.append(missing('unseen_reference_bias','integration',
                    'X1/X2 and anchor discrepancies do not bound a common reference bias; no arbitrary 1e-11 floor'))
                terms.append(missing('reference_parameter_remainders','eos_equilibrium',
                    'finite EOS/center sensitivities and sparse-target transfer are not unseen-remainder bounds'))
                terms.append(term('fit_deterministic_all21','nonlinear_fit',det_fit,'measured_difference'))
                terms.append(conditional('fit_reference_inputs','nonlinear_fit',fit_uncertain,
                    'reference-input-error-enclosure',dependencies=('reference_interpolation','X2_discretization','X2_common_B_alignment',
                    'X2_EOS_sensitivity','X2_center_sensitivity','tangent_stencil' if n.startswith('t_') else 'direct_nonlinear_quantity')))
                terms.append(term('fit_design_floor_observations','nonlinear_fit',fit_design,'measured_difference'))
                terms.append(missing('fit_nonlinear_input_remainder','nonlinear_fit',
                    'design and value-dependent weight perturbations beyond retained finite differences; correlated B0 input remainder'))
                gap=max((K[n]*g['tested_spacing']**2/8 for g in curve['intervals'][:i]),default=0.)
                terms.append(conditional('physical_inter_sample_remainder','between_samples',gap,
                    'global-curvature-and-cut-law') if continuation_checks else
                    missing('physical_inter_sample_remainder','between_samples','frozen continuation sample checks failed'))
                numerator=h['drift_bound']+sum(t['value'] or 0 for t in terms)
                terms.append(zero('absolute_tangent_normalization','normalization','historical absolute tangent units S=1') if n.startswith('t_') else
                    conditional('full_numerator_normalization','normalization',numerator*u0/(1-u0),
                    'B0-reference-normalization',dependencies=tuple(t['name'] for t in terms)))
                if n=='mu_n_profile':
                    terms.append(conditional('equilibrium_profile_definition','production_representation',0.,
                                             'equilibrium-mu-n-profile-identity'))
                m1_new[f'{i}:{n}']=dict(discrepancy=float(D[i]),partial_reference_band=known_band,
                    partial_R=float(D[i])+known_band,share=.23*(h['threshold']-h['drift_bound']),
                    status=m1_decision(float(D[i]),None,.23*(h['threshold']-h['drift_bound']),complete=False,accepted=False),
                    reference_band_complete=False,resolution_ratio=known_band/float(D[i]) if D[i] else None)
            elif n in ('P2_average','species_support'):
                if not np.all(P==0):
                    raise RuntimeError('historical exact identity violated')
                terms=[zero('exact_'+c,c,'original P2 identity' if n=='P2_average' else 'sampled topology identity') for c in FRACTIONS]
                if n=='species_support':
                    terms=[t for t in terms if t['category']!='between_samples']
                    terms.append(missing('support_continuation','between_samples',
                        'fixed threshold ordering and positive B_H continuation have not been accepted as a continuous enclosure'))
            else:
                # R remains a separate incomplete route for metric_structure.
                reconstruction=np.stack([arr(comp['R'][str(j)][n]) for j in range(21)])
                rep=P-reconstruction
                representation=float(np.max(abs(rep[:i+1]-rep[0])))
                terms=[missing('unresolved_'+c,c,'nodewise X2 reference absent; no dense-output accuracy/coverage benchmark')
                       for c in FRACTIONS if c!='production_representation']
                terms.append(term('deterministic_P_minus_R','production_representation',representation,'measured_difference'))
                terms.append(missing('nodewise_R_minus_Q','production_representation',
                    'P-R is not P-Q; metric observer withheld because reference/continuation/authority prerequisites remain unresolved'))
            decision=qualify(h,terms,registry,frozen_historical=oldrows[(i,n)]['historical'],
                             exact=n in ('P2_average','species_support'))
            rows.append(dict(index=i,quantity=n,historical=h,historical_assessment_status=oldrows[(i,n)]['status'],
                headroom=h['threshold']-h['drift_bound'],units=oldrows[(i,n)]['units'],component_count=P.shape[1],
                terms=terms,reference='X2 at achieved B' if n in refs else 'exact identity' if n in ('P2_average','species_support') else
                'R diagnostic only; independent nodewise reference missing',**decision))
    counts={s:sum(r['status']==s for r in rows) for s in ('PASS','FAIL','UNQUALIFIED')}
    summary=dict(schema='compactstar.ba13.evidence-completion.v1',counts=counts,owner_ratified=False,
        certificate_promoted=False,entry_commit='093d121135e35bc50ca2a807470906bcb28ef5c3',
        predeclaration_commit='93c28c3bd85afea8a893a34a8a6dd92ce06c52be',
        historical_outcomes_preserved=True,full_BA13_qualified=False,
        numerical_expansion_stopped=True,stop_reason='reference unseen remainders, continuous continuation and evidence authority unresolved',
        milestones=dict(tangents=tangents['status'],holdouts=holdouts['status'],
            continuous='CONDITIONAL',mu_n_profile='IDENTITY_SUPPORTED_DEFINITION_PENDING',
            metric_structure='NODEWISE_REFERENCE_UNQUALIFIED',qualification='UNQUALIFIED'),
        certificate_sha256=sha(R/'ba13-recovery/run1/fit/certificate.json'),
        fitter_sha256=sha(S/'tests/bnv/compare_candidate.py'),
        acceptance_registry_sha256=sha(S/'docs/validation/ba13-evidence/acceptance-registry.json'))
    def dump(name,obj):
        with (O/name).open('x') as f:
            json.dump(obj,f,indent=1,allow_nan=False,default=lambda v:float(v) if isinstance(v,np.floating) else v.tolist())
            f.write('\n')
    dump('matrix.json',dict(**summary,rows=rows));dump('summary.json',summary)
    dump('comparison.json',details);dump('successor-m1.json',m1_new)
    dump('continuation.json',dict(window=m3['window'],moving_surface=diag['E_surface'],
         groups={n:{k:v for k,v in d.items() if k in ('K','curvature_ratio','continuation_sample_checks','validated_supremum')}
                 for n,d in details.items()},
         physical_response='fixed topology analytic equilibrium response; snapped numerical support/envelope nodes are representation choices',
         unresolved=['uniform curvature supremum','positive B_H away from sampled stencils','continuous surface-cut law'],
         observer=dict(needed_for_full_metric_claim=True,implemented=False,benchmark_run=False,
            reason='reference uncertainty and continuation authority stop; quadratic dense output order is not endpoint GL4 order')))
    cols=['index','quantity','historical_status','historical_utilization','historical_headroom','finite_added_estimates',
          'added_justified_error','remaining_headroom_after_finite_terms','remaining_physical_headroom','status','unresolved_or_conditional']
    with (O/'matrix.tsv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader()
        for r in rows:
            w.writerow(dict(index=r['index'],quantity=r['quantity'],historical_status=r['historical_assessment_status'],
                historical_utilization=r['historical']['utilization'],historical_headroom=r['headroom'],
                **{k:r[k] for k in cols if k in r and k not in ('index','quantity')},
                unresolved_or_conditional='; '.join(r['refusal_reasons'])))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
