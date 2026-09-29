from pathlib import Path
import hashlib,json
q=Path(__file__).resolve().parent;src=Path('/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration')
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def normalize_ranlib_timestamp(path):
 data=bytearray(path.read_bytes());assert data[:8]==b'!<arch>\n';offset=8;members=0
 while offset<len(data):
  header=data[offset:offset+60];assert header[58:]==b'`\n';size=int(header[48:58]);name=bytes(header[:16]).decode().strip()
  if name.startswith('#1/'):
   length=int(name[3:]);name=bytes(data[offset+60:offset+60+length]).rstrip(b'\0').decode()
  if name=='__.SYMDEF':data[offset+16:offset+28]=b'0           '
  offset+=60+size+(size%2);members+=1
 return bytes(data),members
reports={}
for mode in ['Debug','Release']:
 prefix=q/f'install-{mode}'
 assert normalize_ranlib_timestamp(prefix/'lib/libCompactStar.a')==normalize_ranlib_timestamp(q/f'new-{mode}/libCompactStar.a')
 expected={str(p.relative_to(src)):digest(p) for p in (src/'CompactStar').rglob('*') if p.is_file() and p.suffix in ['.h','.hpp'] and p.name not in ['CompactStarConfig.h','CompactStarConfig 2.h']}
 expected['CompactStar/Core/CompactStarConfig.h']=digest(q/f'new-{mode}/generated/include/CompactStar/Core/CompactStarConfig.h')
 actual={str(p.relative_to(prefix/'include')):digest(p) for p in (prefix/'include').rglob('*') if p.is_file()}
 assert actual==expected,(mode,set(actual)^set(expected))
 reports[mode]=dict(build_library_sha256=digest(q/f'new-{mode}/libCompactStar.a'),archive_member_count=normalize_ranlib_timestamp(prefix/'lib/libCompactStar.a')[1],archive_difference='Only the __.SYMDEF timestamp refreshed by CMake install ranlib; every member payload and all other bytes are exact.',library_sha256=digest(prefix/'lib/libCompactStar.a'),installed_header_count=len(actual),installed_header_manifest_sha256=hashlib.sha256(''.join(f'{v}  {k}\n' for k,v in sorted(actual.items())).encode()).hexdigest(),all_public_headers_match=True,historical_headers_installed=False,prefix=str(prefix))
(q/'evidence/final-install-audit.json').write_text(json.dumps(reports,indent=2)+'\n')
print('Fresh installed libraries and complete public-header trees match both final builds')
