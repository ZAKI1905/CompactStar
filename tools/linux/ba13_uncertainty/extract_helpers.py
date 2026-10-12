"""Extract unchanged historical validation helpers, excluding its CLI main."""
import hashlib,sys
from pathlib import Path
src=Path(sys.argv[1]);out=Path(sys.argv[2]);b=src.read_bytes()
assert hashlib.sha256(b).hexdigest()=='30d59d0026902246a2e96778914b6a79f0aaf4a10ac76d519390d13bba314e80'
prefix=b.decode().split('\nint main(int argc,char** argv)')[0]
prefix=prefix.replace('#include "../analysis/phase5b_freegas_validation.cpp"','#include "'+str(src.parent.parent/'analysis/phase5b_freegas_validation.cpp')+'"')
out.write_text(prefix+'\n')
