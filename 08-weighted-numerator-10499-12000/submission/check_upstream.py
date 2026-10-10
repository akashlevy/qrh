import os,json,subprocess,collections,sys
from lhash import ofString, ofBytes
B='/Users/akashlevy/Documents/riemann/formalize_11_12/upstream/.lake/build/lib/lean'
GIT='/Users/akashlevy/Documents/riemann/formalize_11_12/math'
tracked=set(subprocess.check_output(['git','-C',GIT,'ls-files','lean']).decode().split('\n'))
def flat(x,out):
    if isinstance(x,list):
        if len(x)==2 and isinstance(x[0],str) and isinstance(x[1],str): out.append((x[0],x[1])); return
        if len(x)==2 and isinstance(x[0],str): out.append((x[0],None)); flat(x[1],out); return
        for y in x: flat(y,out)
res=collections.Counter(); bad=collections.defaultdict(list); info={}
for dp,dn,fn in os.walk(B):
    for f in fn:
        if not f.endswith('.olean'): continue
        p=os.path.join(dp,f); m=os.path.relpath(p,B)[:-6].replace('/','.')
        tp=p[:-6]+'.trace'
        if not os.path.exists(tp): bad['no_trace'].append(m); continue
        t=json.load(open(tp)); items=[]; flat(t['inputs'],items)
        srcs=[(k,v) for k,v in items if k.endswith('.lean')]
        if len(srcs)!=1: bad['nsrc'].append(m); continue
        sp,sh=srcs[0]
        rp=os.path.normpath(sp)
        exp=os.path.join(GIT,'lean',m.replace('.','/')+'.lean')
        if rp!=exp: bad['srcpath_other'].append((m,rp))
        rel=os.path.relpath(exp,GIT)
        if rel not in tracked: bad['src_not_tracked_at_HEAD'].append(m); 
        else:
            data=subprocess.run(['git','-C',GIT,'show','HEAD:'+rel],capture_output=True).stdout
            if '%016x'%ofString(data.replace(b'\r\n',b'\n'))!=sh: bad['src_hash_mismatch_HEAD'].append(m)
            else: res['src_match_HEAD']+=1
        o=t['outputs']['o']; o=o[0] if isinstance(o,list) else o
        oh=o.split('.')[0]
        if '%016x'%ofBytes(open(p,'rb').read())!=oh: bad['olean_hash_mismatch'].append(m)
        else: res['olean_match']+=1
        hf=p+'.hash'
        if os.path.exists(hf) and open(hf).read().strip()!=oh: bad['hashfile_mismatch'].append(m)
        imps=[k.split(':')[0] for k,v in items if k.endswith(':importAllArts')]
        lean=[k for k,v in items if k.startswith('Lean 4')]
        info[m]={'imports':imps,'lean':lean}
        res['total']+=1
print(dict(res))
for k,v in bad.items(): print(k,len(v),v[:8])
json.dump(info,open('upstream_info.json','w'))
print(collections.Counter(tuple(v['lean']) for v in info.values()))
