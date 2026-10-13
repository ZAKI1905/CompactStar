import sys,json
from pathlib import Path
from equilibrium import analyze
S=Path(__file__).resolve().parents[3];R=S.parents[1];out=R/'ba13-onset-pilot'
case=sys.argv[1];index=int(sys.argv[2]);source=R/('ba13-uncertainty' if len(sys.argv)>3 and sys.argv[3]=='old' else 'ba13-onset-pilot')/'runs'/f'{case}-{index}'
hp=json.loads((out/'constants.json').read_text());d=analyze(source,hp,extra=(case=='r640000'))
(out/'analysis').mkdir(exist_ok=True);(out/'analysis'/f'{case}-{index}.json').write_text(json.dumps(d,indent=2)+'\n')
