import gzip, os
fx='sources/CompactStar/tests/dependency_migration/fixtures/'
base='NStar/Dark_Core/0.8/B_conts/2.01/'
def cmp(rel, label):
    try:
        a=gzip.open(fx+'Debug/'+rel+'.gz','rt').read().splitlines(); b=gzip.open(fx+'Release/'+rel+'.gz','rt').read().splitlines()
    except FileNotFoundError as e: print(label,'missing',e); return
    worst=0; nd=0; ex=None
    for i,(x,y) in enumerate(zip(a,b)):
        if x==y: continue
        tx,ty=x.split('\t'),y.split('\t')
        for j,(u,v) in enumerate(zip(tx,ty)):
            if u!=v:
                try: fu,fv=float(u),float(v); r=abs(fu-fv)/max(abs(fu),abs(fv),1e-300)
                except ValueError: r=0
                nd+=1
                if r>worst: worst=r; ex=(i,j,u.strip(),v.strip())
    print('%-34s rows=%d/%d differing_tokens=%d worst_rel=%.2e example=%s' % (label,len(a),len(b),nd,worst,ex))
print('WITHIN-MAC CONTROL: Mac Debug fixture vs Mac Release fixture')
cmp(base+'BNV_tau/BNV_tau_neutron.tsv','BNV_tau_neutron')
cmp(base+'BNV_tau/BNV_tau_lambda.tsv','BNV_tau_lambda')
cmp(base+'BNV_tau/BNV_tau_sigmam.tsv','BNV_tau_sigmam')
for k in range(10): cmp(base+'2.01_%d_Sequence.tsv'%k,'2.01_%d_Sequence'%k)
cmp('NStar/Dark_Core/0.8/Crit_curve_smooth.tsv','Crit_curve_smooth')
cmp('NStar/Dark_Core/0.8/0.8_19x19_Sequence.tsv','19x19_Sequence')
cmp('EOS/Fermi_Gas_0.8mn.eos','Fermi_Gas_EOS')
