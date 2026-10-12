import gzip, struct, sys, glob, os
def ulps(a,b):
    ia=struct.unpack('<q',struct.pack('<d',a))[0]; ib=struct.unpack('<q',struct.pack('<d',b))[0]; return abs(ia-ib)
def floats(lines):
    out=[]
    for l in lines:
        if l.startswith('#'): continue
        for tok in l.replace(',', ' ').split():
            try: out.append(float(tok))
            except ValueError: out.append(tok)
    return out
def compare(name, mac_lines, lin_lines):
    a, b = floats(mac_lines), floats(lin_lines)
    if len(a) != len(b): print(name, 'TOKEN COUNT DIFFERS', len(a), len(b)); return
    d = [(x,y) for x,y in zip(a,b) if x != y and isinstance(x,float) and isinstance(y,float)]
    nonnum = sum(1 for x,y in zip(a,b) if x != y and not (isinstance(x,float) and isinstance(y,float)))
    if not d: print('%-32s tokens=%d identical (nonnumeric diffs=%d)' % (name, len(a), nonnum)); return
    print('%-32s tokens=%d diff=%d max_rel=%.2e max_ulp=%d nonnumeric_diffs=%d' % (name, len(a), len(d), max(abs(x-y)/max(abs(x),abs(y)) for x,y in d), max(ulps(x,y) for x,y in d), nonnum))
root = 'sources/CompactStar/tests/dependency_migration/fixtures/Debug'
run = glob.glob('toolchain/tmp/csm-fe6_iw1i/run*')[0] if glob.glob('toolchain/tmp/csm-fe6_iw1i/run*') else 'toolchain/tmp/csm-fe6_iw1i/run'
print('linux run root:', run)
for rel in ['EOS/Fermi_Gas_0.8mn.eos', 'NStar/Dark_Core/0.8/B_conts/2.01/BNV_tau/BNV_tau_neutron.tsv', 'NStar/Dark_Core/0.8/B_conts/2.01/BNV_tau/BNV_tau_lambda.tsv', 'NStar/Dark_Core/0.8/B_conts/2.01/BNV_tau/BNV_tau_sigmam.tsv', 'NStar/Dark_Core/0.8/Crit_curve_smooth.tsv', 'NStar/Dark_Core/0.8/0.8_19x19_Sequence.tsv']:
    mac = os.path.join(root, rel + '.gz'); lin = os.path.join(run, rel)
    if not os.path.exists(mac) or not os.path.exists(lin):
        print(rel, 'missing:', os.path.exists(mac), os.path.exists(lin)); continue
    compare(rel.split('/')[-1], gzip.open(mac, 'rt').read().splitlines(), open(lin).read().splitlines())
