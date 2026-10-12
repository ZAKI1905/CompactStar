import gzip, os, sys
root='sources/CompactStar/tests/dependency_migration/fixtures/Debug'; run='toolchain/tmp/csm-fe6_iw1i/run'
def rows(txt): return [l for l in txt.splitlines()]
for rel, maxshow in [('NStar/Dark_Core/0.8/B_conts/2.01/BNV_tau/BNV_tau_neutron.tsv', 60), ('EOS/Fermi_Gas_0.8mn.eos', 14), ('NStar/Dark_Core/0.8/Crit_curve_smooth.tsv', 8), ('NStar/Dark_Core/0.8/0.8_19x19_Sequence.tsv', 6)]:
    mac = rows(gzip.open(os.path.join(root, rel+'.gz'),'rt').read()); lin = rows(open(os.path.join(run, rel)).read())
    print('=====', rel, 'rows', len(mac), len(lin)); print('header:', mac[0][:200])
    shown=0
    for i,(a,b) in enumerate(zip(mac,lin)):
        if a==b: continue
        ta, tb = a.split('\t'), b.split('\t')
        if len(ta)!=len(tb): print('row',i,'TOKEN COUNT DIFFERS'); continue
        for j,(x,y) in enumerate(zip(ta,tb)):
            if x!=y:
                try: fx,fy=float(x),float(y); rel_=abs(fx-fy)/max(abs(fx),abs(fy),1e-300)
                except ValueError: rel_=float('nan')
                if shown<maxshow: print('  row %4d col %2d mac=%-24s linux=%-24s rel=%.2e' % (i,j,x,y,rel_)); shown+=1
