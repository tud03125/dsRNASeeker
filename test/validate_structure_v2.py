"""Real-tool regression; same v2 function, exact ordered pair IDs and output schema."""
import argparse, json, os, sys, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from modules import energetics as serial, structure_execute as opt
from modules.intarna_checked import preflight_intarna

def main():
    p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--out',required=True);p.add_argument('--intarna',required=True);p.add_argument('--n',type=int,default=32);a=p.parse_args()
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    ps=opt.pairs_in(a.fasta); ids=sorted(ps,key=lambda p:(len(ps[p]['A'])+len(ps[p]['B']),p));sel=[ids[i] for i in np.linspace(0,len(ids)-1,min(a.n,len(ids)),dtype=int)]
    args=SimpleNamespace(rnacofold_exe='RNAcofold',rnafold_exe='RNAfold',intarna_exe=a.intarna,threads=1)
    report={'intarna':preflight_intarna(a.intarna),'sample_n':len(sel),'comparisons':{}}
    for cname,mod,cond in [('serial',serial,'CASE'),('parallel',opt,'CASE'),('parallel',opt,'CONTROL')]:
        d=out/cname/cond;d.mkdir(parents=True);fa=d/f'duplex_arms.inverted.{cond}.clean.fa';opt.write_fa(fa,[(p,ps[p]) for p in sel]);t=time.monotonic()
        _,cf=mod.run_rnacofold(args,fa,d,'inverted',cond);mod.run_ddg(args,fa,cf,d,'inverted',cond);mod.run_interface_bpp(fa,d,'inverted',cond)
        report[cname+'_'+cond+'_seconds']=time.monotonic()-t
    for stage in ['cofold','ddg','bpp']:
        def get(c,cond):
            _,f,_=opt.paths(out/c/cond,'inverted',cond,stage);return pd.read_csv(f,sep='\t').set_index('pair_id').loc[sel]
        x=get('serial','CASE');result={}
        for c,cond in [('parallel','CASE'),('parallel','CONTROL')]:
            y=get(c,cond);assert x.columns.equals(y.columns);assert x.index.equals(y.index)
            np.testing.assert_array_equal(x.isna(),y.isna()); assert np.isfinite(x.to_numpy(dtype=float)).all() and np.isfinite(y.to_numpy(dtype=float)).all()
            np.testing.assert_allclose(x,y,atol=1e-10,rtol=0)
            for col in ['lenA','lenB','len_total','interface_bpp_n','interface_bpp_n_ge_1e5']:
                if col in x:np.testing.assert_array_equal(x[col],y[col])
            result[c+'_'+cond]=float(np.max(np.abs(x.to_numpy()-y.to_numpy())))
        _,_,m=opt.paths(out/'parallel'/'CONTROL','inverted','CONTROL',stage)
        result['control_cache']=json.loads(m.read_text());assert result['control_cache']['reused']+result['control_cache']['computed']==len(sel)
        if stage in ('cofold','bpp'): assert result['control_cache']['reused']==len(sel)
        # ddG may conservatively recompute finite rows if embedded MFE floats differ; compare outputs above.
        report['comparisons'][stage]=result
    (out/'STRUCTURE_V2_PASS.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS: production v2 structural smoke',flush=True)
if __name__=='__main__':main()
