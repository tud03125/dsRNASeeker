"""Real ViennaRNA regression; execute from the candidate repository root."""
import argparse, json, os, sys, tempfile, time
from pathlib import Path
from types import SimpleNamespace
import pandas as pd
import numpy as np
from modules import energetics as e

def main():
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,required=True);a.add_argument('--n',type=int,default=50);a=a.parse_args();root=a.out.resolve();root.mkdir(parents=True,exist_ok=False)
    pairs={'p1':('AUGCUGCAUGGCAUUGCUAGCU','AGCUAGCAAUGCCAUGCAGCAU'), 'p2':('AACCGGUUAACCGGUUACGU','ACGUAACCGGUUAACCGGUU'), 'p3':('GCGCGCGCGCGC','GCGCGCGCGCGC'), 'p4':('AUUGCGUACACGUGUAACGUACGU','ACGUACGUUACACGUGUACGCAAU')}
    args=SimpleNamespace(null_n=a.n,null_seed=1,rnafold_exe='RNAfold',rnacofold_exe='RNAcofold',threads=4)
    def inputs(d, ps, shift=0):
        d.mkdir();fa=d/'arms.fa';obs=d/'ddg.tsv';fa.write_text(''.join(f'>{p}|A\n{x}\n>{p}|B\n{y}\n' for p,(x,y) in ps.items()));pd.DataFrame([{'pair_id':p,'ddG_interaction_kcalmol':-20.0+shift} for p in ps]).to_csv(obs,sep='\t',index=False);return fa,obs
    def run(name,fn,ps=pairs,shift=0):
        d=root/name;f,o=inputs(d,ps,shift);t=time.monotonic();r=fn(args,f,o,d,'inverted','case');return pd.read_csv(r,sep='\t').set_index('pair_id').sort_index(),time.monotonic()-t,r
    def equal(x,y):pd.testing.assert_frame_equal(x,y,check_exact=False,atol=1e-10,rtol=0)
    def meta(p):return json.loads(p.with_suffix('.tsv.provenance.json').read_text())
    os.environ['DSRNA_NULLZ_CACHE']=str(root/'cache');os.environ['DSRNA_NULLZ_WORKERS']='4'
    old,t0,_=run('serial',e._run_null_z_legacy);new,t1,f=run('parallel',e.run_null_z);equal(old,new);assert meta(f)['computed']==len(pairs)
    cached,t2,f=run('cached',e.run_null_z,dict(reversed(list(pairs.items()))));equal(old,cached);assert meta(f)['reused']==len(pairs)
    # Remove one checkpoint: restart must recalculate just that pair.
    next((root/'cache').rglob('*.json')).unlink();x,_,f=run('resumed',e.run_null_z);equal(old,x);assert meta(f)['computed']==1
    args.null_seed=2;x,_,f=run('seed_changed',e.run_null_z);assert meta(f)['computed']==len(pairs);y,_,_=run('seed_serial',e._run_null_z_legacy);equal(x,y)
    args.null_seed=1;args.null_n=a.n+1;x,_,f=run('n_changed',e.run_null_z);assert meta(f)['computed']==len(pairs);y,_,_=run('n_serial',e._run_null_z_legacy);equal(x,y)
    args.null_n=a.n;x,_,f=run('obs_changed',e.run_null_z,shift=1);assert meta(f)['computed']==len(pairs);y,_,_=run('obs_serial',e._run_null_z_legacy,shift=1);equal(x,y)
    altered=dict(pairs);altered['p1']=('A'+pairs['p1'][0],pairs['p1'][1]);x,_,f=run('seq_changed',e.run_null_z,altered);assert meta(f)['computed']==1;y,_,_=run('seq_serial',e._run_null_z_legacy,altered);equal(x,y)
    report={'pairs':len(pairs),'n':a.n,'serial_seconds':t0,'parallel_seconds':t1,'cached_seconds':t2,'checks':['serial/parallel equivalence','cache reorder','partial checkpoint resume','seed invalidation','N invalidation','observed ddG invalidation','sequence invalidation'],'tolerance':1e-10}
    (root/'NULLZ_PASS.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS: Null-Z equivalence and cache invalidation',report)
if __name__=='__main__':main()
