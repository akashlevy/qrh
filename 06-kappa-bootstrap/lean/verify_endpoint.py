"""Check the public overlay endpoint and report its actual kernel axioms."""
from pathlib import Path
import json,os,subprocess,time
out=Path(__file__).resolve().parent
conf=json.loads((out/'lean_environment.json').read_text())
env=os.environ.copy();env['LEAN_PATH']=conf['LEAN_PATH']
log=out/'logs/endpoint-axioms.log'
start=time.time()
with log.open('w') as stream:
 result=subprocess.run([conf['lean'],'-DautoImplicit=false','-M','8192','-j','2','-R',str(out),'-o',str(out/'KappaFormal.olean'),str(out/'KappaFormal.lean')],cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1200)
text=log.read_text()
print(text)
status={'exit_code':result.returncode,'elapsed_seconds':time.time()-start,'has_sorryAx':'sorryAx' in text,'endpoint_log':str(log)}
(out/'endpoint_verification.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status))
raise SystemExit(result.returncode)
