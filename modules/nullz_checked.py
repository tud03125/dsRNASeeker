"""Parallel/resumable execution of the unchanged v2 Null-Z computation.
Only complete, fingerprinted per-pair results are reusable. No shared temp files.
"""
from __future__ import annotations
import hashlib, json, math, os, shutil, subprocess, sys, tempfile, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import pandas as pd
from . import energetics as legacy

COLS = ['pair_id','ddG_mu_null','ddG_sd_null','ddG_Z','null_n_requested','null_n_effective','null_shuffle_exact_dinuc','null_shuffle_method']
METHOD='randomized_eulerian_exact_dinucleotide'
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()
def digest(x):
    return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def write_atomic(path, data):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f:
            f.write(data);f.flush();os.fsync(f.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
def executable_info(exe):
    p=shutil.which(str(exe))
    if not p: raise RuntimeError(f'Folding executable not found: {exe}')
    p=Path(p).resolve(); ver=subprocess.run([str(p),'--version'],capture_output=True,text=True,check=True,timeout=30)
    # Include dynamically loaded libraries where Linux ldd exposes them.
    libs={}
    if shutil.which('ldd'):
        q=subprocess.run(['ldd',str(p)],capture_output=True,text=True,timeout=30)
        if 'not found' in q.stdout: raise RuntimeError(f'Unresolved libraries: {p}: {q.stdout}')
        for ln in q.stdout.splitlines():
            for tok in ln.split():
                if tok.startswith('/') and Path(tok).is_file():
                    lib=Path(tok).resolve();libs[str(lib)]=sha(lib)
    return {'path':str(p),'sha256':sha(p),'version':ver.stdout.strip(),'libraries':libs}
def fingerprint(args):
    return {'schema':1,'legacy_source':sha(legacy.__file__),'wrapper_source':sha(__file__),
            'python':sys.version,'pandas':pd.__version__, 'numpy':legacy.np.__version__,
            'fold':executable_info(args.rnafold_exe),'cofold':executable_info(args.rnacofold_exe),
            'environment':{k:os.environ.get(k) for k in ('CONDA_PREFIX','LD_LIBRARY_PATH','VRNA_PARAMS_PATH','OMP_NUM_THREADS','LC_ALL','LANG')},
            'n':int(args.null_n),'seed':int(args.null_seed),'fold_options':['--noPS']}
def read_pairs(path):
    pairs={}; seen=set(); key=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            pid,arm=line[1:].strip().rsplit('|',1)
            if not pid or arm not in ('A','B') or (pid,arm) in seen: raise ValueError(f'Malformed or duplicate FASTA record: {line}')
            key=(pid,arm);seen.add(key);pairs.setdefault(pid,{})[arm]=''
        elif line.strip():
            if key is None: raise ValueError('Sequence before header')
            pairs[key[0]][key[1]]+=line.strip()
    if any(set(ab)!= {'A','B'} or not all(ab.values()) for ab in pairs.values()): raise ValueError('Incomplete FASTA arms')
    return pairs

def validate_row(row,pid,n):
    if set(row)!=set(COLS) or row['pair_id']!=pid or row['null_n_requested']!=n or row['null_n_effective']!=n or row['null_shuffle_exact_dinuc'] is not True or row['null_shuffle_method']!=METHOD:
        raise ValueError(f'Incomplete or incompatible Null-Z result: {pid}')
    mu,sd,z=(row[k] for k in ('ddG_mu_null','ddG_sd_null','ddG_Z'))
    if not isinstance(mu,(float,int)) or not math.isfinite(mu): raise ValueError(f'Invalid null mean: {pid}')
    if n==1:
        if sd is not None or z is not None: raise ValueError(f'Invalid single-replicate null: {pid}')
    elif not isinstance(sd,(float,int)) or not math.isfinite(sd) or sd<0: raise ValueError(f'Invalid null SD: {pid}')
    elif sd==0:
        if z is not None: raise ValueError(f'Zero variance requires missing Z: {pid}')
    elif z is None or not math.isfinite(z): raise ValueError(f'Invalid Z: {pid}')
    return row

def read_cache(path,key,pid,n):
    if not path.is_file():return None
    try:
        d=json.loads(path.read_text());r=d['row']
        if d['key']!=key or d['row_sha256']!=digest(r):return None
        return validate_row(r,pid,n)
    except (KeyError,ValueError,TypeError,OSError): return None

def compute(args,pid,ab,obs,cache,key,tag,cond):
    n=int(args.null_n)
    hit=read_cache(cache,key,pid,n)
    if hit is not None:return hit,True
    with tempfile.TemporaryDirectory(prefix='dsrna_nullz_',dir=os.environ.get('TMPDIR')) as t:
        t=Path(t);fa=t/'pair.fa'; ddg=t/'obs.tsv'
        fa.write_text(f'>{pid}|A\n{ab["A"]}\n>{pid}|B\n{ab["B"]}\n')
        pd.DataFrame([{'pair_id':pid,'ddG_interaction_kcalmol':obs}]).to_csv(ddg,sep='\t',index=False)
        out=legacy._run_null_z_legacy(args,fa,ddg,t,tag,cond)
        df=pd.read_csv(out,sep='\t',dtype={'pair_id':str},float_precision='round_trip')
        if len(df)!=1 or list(df.columns)!=COLS:raise ValueError(f'Bad legacy output: {pid}')
        v=df.iloc[0];row={'pair_id':pid,'null_n_requested':int(v.null_n_requested),'null_n_effective':int(v.null_n_effective),'null_shuffle_exact_dinuc':bool(v.null_shuffle_exact_dinuc),'null_shuffle_method':str(v.null_shuffle_method)}
        for k in ('ddG_mu_null','ddG_sd_null','ddG_Z'):
            row[k]=None if pd.isna(v[k]) else float(v[k])
        validate_row(row,pid,n)
        write_atomic(cache,json.dumps({'key':key,'row':row,'row_sha256':digest(row)},sort_keys=True,allow_nan=False)+'\n')
        return row,False

def run_nullz_checked(args,clean_fa,ddg_tsv,outdir,tag,condition):
    outdir=Path(outdir).resolve();outdir.mkdir(parents=True,exist_ok=True)
    n=int(args.null_n)
    if n<1:raise ValueError('--null-n must be >=1')
    pairs=read_pairs(clean_fa);obs=pd.read_csv(ddg_tsv,sep='\t',dtype={'pair_id':str},float_precision='round_trip').set_index('pair_id',verify_integrity=True)
    if set(obs.index)!=set(pairs):raise ValueError('Null-Z FASTA/ddG pair IDs differ; cannot certify inputs')
    vals=pd.to_numeric(obs.ddG_interaction_kcalmol,errors='raise')
    if not all(math.isfinite(v) for v in vals):raise ValueError('Nonfinite observed ddG; inspect upstream calculation')
    fp=fingerprint(args);fp_hash=digest(fp)
    cache=Path(os.environ.get('DSRNA_NULLZ_CACHE',str(outdir.parent/'.dsrna_nullz_cache'))).resolve();cache.mkdir(parents=True,exist_ok=True)
    # Absolute executable paths keep worker subprocesses independent of cwd.
    from types import SimpleNamespace
    local=SimpleNamespace(null_n=n,null_seed=int(args.null_seed),rnafold_exe=fp['fold']['path'],rnacofold_exe=fp['cofold']['path'])
    cpus=int(os.environ.get('SLURM_CPUS_PER_TASK',getattr(args,'threads',1) or 1));workers=max(1,min(cpus,int(os.environ.get('DSRNA_NULLZ_WORKERS',str(cpus)))))
    out=outdir/f'duplex_pairs.{tag}.{condition}.nullZ.tsv'; manifest=out.with_suffix('.tsv.provenance.json')
    # Old aggregate output is never trusted; preserve it until successful replacement.
    start=time.monotonic();rows={};reused=0
    tasks=[]
    for pid in obs.index:
        ab=pairs[pid];v=float(vals.loc[pid]);key=digest([fp_hash,pid,ab['A'],ab['B'],v]);p=cache/key[:2]/(key+'.json')
        hit=read_cache(p,key,pid,n)
        if hit is None:tasks.append((pid,ab,v,p,key))
        else:rows[pid]=hit;reused+=1
    print(f'[NULL-Z] {condition}: total={len(pairs)} reused={reused} compute={len(tasks)} workers={workers}',flush=True)
    tasks.sort(key=lambda x:len(x[1]['A'])+len(x[1]['B']),reverse=True)
    # Bounded submissions avoid millions of queued futures. Every pair is an atomic checkpoint.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        it=iter(tasks);pending={}
        def enqueue():
            try:pid,ab,v,p,key=next(it)
            except StopIteration:return False
            fut=pool.submit(compute,local,pid,ab,v,p,key,tag,condition);pending[fut]=pid;return True
        for _ in range(workers*2):
            if not enqueue():break
        done=0
        while pending:
            finished,_=wait(pending,return_when=FIRST_COMPLETED)
            for f in finished:
                pid=pending.pop(f);r,hit=f.result();rows[pid]=r;reused+=int(hit);done+=1
                if done%16==0 or done==len(tasks):print(f'[NULL-Z] {done}/{len(tasks)} pending pairs finished; {time.monotonic()-start:.1f}s',flush=True)
                enqueue()
    if set(rows)!=set(obs.index):raise RuntimeError('Null-Z output IDs differ')
    df=pd.DataFrame([rows[p] for p in obs.index],columns=COLS)
    write_atomic(out,df.to_csv(sep='\t',index=False))
    meta={'fingerprint':fp,'fasta_sha256':sha(clean_fa),'ddg_sha256':sha(ddg_tsv),'output_sha256':sha(out),'total':len(rows),'reused':reused,'computed':len(rows)-reused,'elapsed_seconds':time.monotonic()-start,'finite_z':int(df.ddG_Z.notna().sum())}
    write_atomic(manifest,json.dumps(meta,indent=2,sort_keys=True)+'\n');return out
