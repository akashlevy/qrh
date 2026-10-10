"""Build only changed descendants, reusing audited original cached dependencies."""
from pathlib import Path
import json,os,re,subprocess,sys,time
BASE=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
SRC=BASE/'formalize_11_12/math/lean'
UP=BASE/'formalize_11_12/upstream'
LEAN=Path('/Users/akashlevy/.elan/toolchains/leanprover--lean4---v4.34.1/bin/lean')
LAKE=Path('/Users/akashlevy/.elan/bin/lake')
prefix='OAI.NumberTheory.DirichletL.'
manifest=json.loads((OUT/'dependency_manifest.json').read_text())
affected=set(manifest['affected'])
def newmod(m):return 'KappaFormal.'+m.removeprefix(prefix)
def newpath(m):return OUT/Path(newmod(m).replace('.','/')).with_suffix('.lean')
imports={}
for m in sorted(affected):
 original=(SRC/Path(m.replace('.','/')).with_suffix('.lean')).read_text()
 imports[m]=[x for line in original.splitlines() if line.startswith('import ') for x in line.split()[1:]]
 transformed=re.sub(r'^import\s+(\S+)',lambda mt:'import '+(newmod(mt[1]) if mt[1] in affected else mt[1]),original,flags=re.M)
 # Weaken only the explicit kappa lower bounds, not unrelated appearances of 3/4.
 transformed=re.sub(r'3\s*/\s*4(?=\s*≤\s*(?:κ|kappa))','13/18',transformed)
 transformed=re.sub(r'3\s*/\s*4(?=\s*:\s*ℝ\)\s*≤\s*(?:κ|kappa))','13/18',transformed)
 if m.endswith('ReferenceLowBranchGeometry'):
  for old,new in [('16*M/21','23*M/30'),('13*M/14','14*M/15'),('M/21','M/20'),('7/2*ell','10/3*ell'),('M/14','M/15')]:
   transformed=transformed.replace(old,new)
 if m.endswith('LiveClippingDefect'):
  for old,new in [('15*M/28-9*defect/7','8*M/15-13*defect/10'),('16*M/21','23*M/30'),('13*M/14','14*M/15'),('M/21+2*defect/7','M/20+3*defect/10'),('7/2*ell','10/3*ell'),('4*defect/7','3*defect/5'),('11*defect/7','8*defect/5'),('4*e/7','3*e/5'),('11*e/7','8*e/5'),('M/14','M/15')]:
   transformed=transformed.replace(old,new)
 # Propagate reserve assumptions through reference/reflection wrappers.
 for old,new in [('xi+11*defect/7','xi+8*defect/5'),('xi+11*e/7','xi+8*e/5'),('xi+11*(rho/100)/7','xi+8*(rho/100)/5'),('xi≤s.width/14','xi≤s.width/15'),('≤s.width/14','≤s.width/15')]:
  transformed=transformed.replace(old,new)
 dest=newpath(m);dest.parent.mkdir(parents=True,exist_ok=True)
 if not dest.exists():dest.write_text(transformed)
order=[];seen=set()
def visit(m):
 if m in seen:return
 seen.add(m)
 for i in imports[m]:
  if i in affected:visit(i)
 order.append(m)
visit(manifest['target'])
(OUT/'compile_order.json').write_text(json.dumps(order,indent=2))
result=subprocess.run([str(LAKE),'env','printenv','LEAN_PATH'],cwd=UP,text=True,capture_output=True,check=True)
base_path=result.stdout.strip().splitlines()[-1]
env=os.environ.copy();env['LEAN_PATH']=str(OUT)+':'+base_path
(OUT/'lean_environment.json').write_text(json.dumps({'lean':str(LEAN),'LEAN_PATH':env['LEAN_PATH']},indent=2))
logs=OUT/'logs';logs.mkdir(exist_ok=True)
for index,m in enumerate(order,1):
 source=newpath(m);ole=source.with_suffix('.olean');log=logs/(newmod(m)+'.log')
 dependencies=[newpath(i).with_suffix('.olean') for i in imports[m] if i in affected]
 if ole.exists() and ole.stat().st_mtime>source.stat().st_mtime and all(ole.stat().st_mtime>d.stat().st_mtime for d in dependencies):
  print(f'[{index}/{len(order)}] cached {newmod(m)}',flush=True);continue
 start=time.time();print(f'[{index}/{len(order)}] compiling {newmod(m)}',flush=True)
 with log.open('w') as out:
  completed=subprocess.run([str(LEAN),'-DautoImplicit=false','-M','8192','-j','2','-R',str(OUT),'-o',str(ole),str(source)],cwd=OUT,env=env,text=True,stdout=out,stderr=subprocess.STDOUT,timeout=1200)
 print(f'  exit={completed.returncode}, {time.time()-start:.1f}s',flush=True)
 if completed.returncode:
  print(log.read_text(),flush=True);sys.exit(completed.returncode)
print('ALL OVERLAY MODULES COMPILED',flush=True)
