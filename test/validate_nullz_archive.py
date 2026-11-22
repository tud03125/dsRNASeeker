"""Compare unchanged serial and optimized null model on actual archived J2 arm sequences."""
import argparse, json, os, time
from pathlib import Path
from types import SimpleNamespace
import pandas as pd
from modules import energetics as e
from modules.nullz_checked import read_pairs

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fasta',type=Path,required=True);ap.add_argument('--ddg',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();d=a.out.resolve();d.mkdir(parents=True,exist_ok=False)
    pairs=read_pairs(a.fasta);obs=pd.read_csv(a.ddg,sep='\t').set_index('pair_id',verify_integrity=True)
    elig=sorted((p for p,ab in pairs.items() if p in obs.index and 0<len(ab['A']) and 0<len(ab['B']) and len(ab['A'])+len(ab['B'])<=600 and pd.notna(obs.loc[p,'ddG_interaction_kcalmol'])),key=lambda p:(len(pairs[p]['A'])+len(pairs[p]['B']),p))
    if len(elig)<4:raise ValueError('Fewer than four archived eligible pairs')
    k=min(8,len(elig));ids=[elig[round(i*(len(elig)-1)/(k-1))] for i in range(k)]
    fa=d/'arms.fa';ddg=d/'obs.tsv';fa.write_text(''.join(f'>{p}|{arm}\n{pairs[p][arm]}\n' for p in ids for arm in ('A','B')));obs.loc[ids].reset_index().to_csv(ddg,sep='\t',index=False)
    args=SimpleNamespace(null_n=50,null_seed=1,rnafold_exe='RNAfold',rnacofold_exe='RNAcofold',threads=4)
    os.environ['DSRNA_NULLZ_CACHE']=str(d/'cache');os.environ['DSRNA_NULLZ_WORKERS']='4'
    def run(name,fn):
        x=d/name;x.mkdir();t=time.monotonic();f=fn(args,fa,ddg,x,'inverted','case');return pd.read_csv(f,sep='\t').set_index('pair_id'),time.monotonic()-t,f
    x,s,_=run('legacy',e._run_null_z_legacy);y,p,_=run('parallel',e.run_null_z);z,c,f=run('cached',e.run_null_z)
    pd.testing.assert_frame_equal(x,y,check_exact=False,atol=1e-10,rtol=0);pd.testing.assert_frame_equal(x,z,check_exact=False,atol=1e-10,rtol=0)
    assert x.null_n_effective.eq(50).all();assert json.loads(f.with_suffix('.tsv.provenance.json').read_text())['reused']==k
    report={'pair_ids':ids,'pair_n':k,'null_n':50,'seed':1,'max_total_arm_length':600,'serial_seconds':s,'parallel_seconds':p,'cached_seconds':c,'tolerance':1e-10,'source_fasta':str(a.fasta),'source_ddg':str(a.ddg)}
    (d/'NULLZ_ARCHIVE_PASS.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS: actual archived sequence Null-Z execution equivalence',report)
if __name__=='__main__':main()
