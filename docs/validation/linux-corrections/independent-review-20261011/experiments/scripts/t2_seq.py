import gzip, os
root='sources/CompactStar/tests/dependency_migration/fixtures/Debug'; run='toolchain/tmp/csm-fe6_iw1i/run'
base='NStar/Dark_Core/0.8/B_conts/2.01/'
print('contour point counts (rows) Mac vs Linux:')
for k in range(10):
    rel=base+'2.01_%d_Sequence.tsv'%k
    m=gzip.open(os.path.join(root,rel+'.gz'),'rt').read().splitlines(); l=open(os.path.join(run,rel)).read().splitlines()
    nd=sum(1 for a,b in zip(m,l) if a!=b)
    print('  k=%d rows mac=%d linux=%d differing_rows=%d' % (k,len(m),len(l),nd))
for rel in [base+'2.01_5_Sequence.tsv', base+'BNV_rates_2.01_5.tsv', base+'2.01_0_Sequence.tsv']:
    m=gzip.open(os.path.join(root,rel+'.gz'),'rt').read().splitlines(); l=open(os.path.join(run,rel)).read().splitlines()
    print('=====',rel.split('/')[-1],'rows',len(m),len(l)); print('header:',m[0][:230])
    shown=0; worst={}
    for i,(a,b) in enumerate(zip(m,l)):
        if a==b: continue
        ta,tb=a.split('\t'),b.split('\t')
        if len(ta)!=len(tb): print(' row',i,'token count differs'); continue
        for j,(x,y) in enumerate(zip(ta,tb)):
            if x!=y:
                try: fx,fy=float(x),float(y); r=abs(fx-fy)/max(abs(fx),abs(fy),1e-300)
                except ValueError: r=float('nan')
                worst[j]=max(worst.get(j,0),r)
                if shown<10: print('  row %3d col %2d mac=%s linux=%s rel=%.2e'%(i,j,x.strip(),y.strip(),r)); shown+=1
    print('  worst rel per column:', {j:'%.1e'%v for j,v in sorted(worst.items())})
