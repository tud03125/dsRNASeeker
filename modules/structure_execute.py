"""Opt-in production v2 execution wrapper. Thermodynamic calculations delegated to the current serial module.
Cache only complete finite rows with exact pair ID and ordered sequences; results in input order.
"""
from __future__ import annotations
import csv, hashlib, json, math, multiprocessing as mp, os, shutil, subprocess, sys, tempfile, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from functools import lru_cache
from . import energetics as serial
from .energetics import prepare_duplex_inputs, run_null_z

SUFFIX = {'cofold':'cofold_mfe','ddg':'ddg','bpp':'interface_bpp','intarna':'IntaRNA'}
COLUMNS = {
 'cofold':['pair_id','RNAcofold_MFE_kcalmol','lenA','lenB','len_total','MFE_norm_kcalpermkb'],
 'ddg':['pair_id','RNAcofold_MFE_kcalmol','lenA','lenB','len_total','MFE_norm_kcalpermkb','RNAfold_A_MFE_kcalmol','RNAfold_B_MFE_kcalmol','ddG_interaction_kcalmol','ddG_norm_kcalpermkb'],
 'bpp':['pair_id','interface_bpp_sum','interface_bpp_max','interface_bpp_n','interface_bpp_n_ge_1e5','interface_bpp_expected_fraction_shorter','interface_bpp_mean_arm_fraction'],
 'intarna':['pair_id','E']}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def pairs_in(fa):
    # Strict parsing avoids silently accepting duplicate or truncated FASTA records.
    pairs={}; key=None; seen=set()
    for line in Path(fa).read_text().splitlines():
        if line.startswith('>'):
            pid,arm=line[1:].strip().rsplit('|',1)
            if not pid or arm not in ('A','B') or (pid,arm) in seen: raise ValueError(f'Invalid/duplicate FASTA record {line}')
            seen.add((pid,arm)); key=(pid,arm); pairs.setdefault(pid,{})[arm]=''
        elif line.strip():
            if key is None: raise ValueError('Sequence before header')
            pairs[key[0]][key[1]]+=line.strip()
    if not pairs or any(set(v)!= {'A','B'} or not v['A'] or not v['B'] for v in pairs.values()):
        raise ValueError(f'Empty or incomplete arms: {fa}; use legacy or inspect inputs')
    return pairs

def write_fa(path, rows):
    with Path(path).open('w') as f:
        for pid,ab in rows:
            for a in ('A','B'): f.write(f'>{pid}|{a}\n{ab[a]}\n')

def read_table(path, stage, ids):
    with Path(path).open(newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        if r.fieldnames!=COLUMNS[stage]: raise ValueError(f'Schema mismatch: {path}: {r.fieldnames}')
        d={}
        for row in r:
            pid=row['pair_id']
            if pid in d or None in row or any(v is None for v in row.values()): raise ValueError(f'Duplicate/malformed row in {path}: {pid}')
            for k,v in row.items():
                if k!='pair_id' and v not in ('','NA','NaN','nan'):
                    float(v)
            d[pid]=row
        if set(d)!=set(ids): raise ValueError(f'Incomplete or mismatched pair IDs: {path}; expected {len(ids)}, observed {len(d)}')
        return d

def write_table(path, stage, rows):
    path=Path(path); tmp=path.with_name(path.name+'.tmp')
    with tmp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=COLUMNS[stage],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    tmp.replace(path)

def write_json(path, d):
    path=Path(path);tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');tmp.replace(path)

def paths(d,tag,condition,stage):
    d=Path(d); fa=d/f'duplex_arms.{tag}.{condition}.clean.fa'; table=d/f'duplex_pairs.{tag}.{condition}.{SUFFIX[stage]}.tsv'
    return fa,table,table.with_suffix(table.suffix+'.cache.json')

@lru_cache(maxsize=8)
def file_info(exe):
    x=shutil.which(exe)
    if not x: raise RuntimeError(f'Executable not found: {exe}')
    x=Path(x).resolve(); return {'path':str(x),'sha256':sha(x)}

def fingerprint(stage,args):
    import RNA
    # Record the Python binding and native extension when exposed by the package.
    native=getattr(getattr(RNA,'_RNA',None),'__file__',None)
    fp={'stage':stage,'legacy_sha256':sha(serial.__file__),'wrapper_sha256':sha(__file__),
        'numpy':serial.np.__version__, 'pandas':serial.pd.__version__, 'python':sys.version,'python_path':str(Path(sys.executable).resolve()),'RNA_version':str(getattr(RNA,'__version__','unknown')),
        'RNA_binding':sha(RNA.__file__),'RNA_native':sha(native) if native and Path(native).is_file() else None,
        'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),'OPENBLAS_NUM_THREADS':os.environ.get('OPENBLAS_NUM_THREADS'),
        'MKL_NUM_THREADS':os.environ.get('MKL_NUM_THREADS')}
    if stage in ('cofold','ddg'): fp['RNAcofold']=file_info(args.rnacofold_exe)
    if stage=='ddg': fp['RNAfold']=file_info(args.rnafold_exe)
    if stage=='intarna': fp['IntaRNA']=file_info(args.intarna_exe)
    return fp

def finite_row(row):
    try:return all(math.isfinite(float(v)) for k,v in row.items() if k!='pair_id')
    except (ValueError,TypeError):return False

def cache_load(d,tag,condition,stage,fp):
    fa,table,manifest=paths(d,tag,condition,stage)
    if not all(p.is_file() for p in (fa,table,manifest)): return None
    meta=json.loads(manifest.read_text())
    if meta['fingerprint']!=fp: return None
    if meta['fasta_sha256']!=sha(fa) or meta['table_sha256']!=sha(table): raise ValueError(f'Cache content changed: {table}')
    pairs=pairs_in(fa); rows=read_table(table,stage,pairs)
    return pairs,rows

def mark_cache(d,tag,condition,stage,fp,extra=None):
    fa,table,manifest=paths(d,tag,condition,stage)
    write_json(manifest,{'fingerprint':fp,'fasta_sha256':sha(fa),'table_sha256':sha(table),**(extra or {})})

def sources(outdir,tag,condition,stage,fp):
    # Current run's other condition (case completes before control in historical workflow).
    for d in sorted(outdir.parent.iterdir()):
        if d.is_dir() and d!=outdir and (d/f'duplex_arms.{tag}.{d.name}.clean.fa').is_file():
            yield d,d.name
    # Explicitly supplied prior production run directories; manifests verified for each table.
    for root in os.environ.get('DSRNA_STRUCT_CACHE_RUNS','').split(os.pathsep):
        if root:
            for d in sorted(Path(root).resolve().glob('*/')):
                if d.resolve()!=outdir.resolve(): yield d,d.name

def work(task):
    stage,args,items,cfrows,tmp,tag,cond=task
    tmp=Path(tmp);tmp.mkdir(parents=True,exist_ok=False)
    fa,table,_=paths(tmp,tag,cond,stage);write_fa(fa,items)
    cwd=Path.cwd()
    try:
        os.chdir(tmp)
        if stage=='cofold': _,ret=serial.run_rnacofold(args,fa,tmp,tag,cond)
        elif stage=='ddg':
            _,ct,_=paths(tmp,tag,cond,'cofold');write_table(ct,'cofold',cfrows)
            ret=serial.run_ddg(args,fa,ct,tmp,tag,cond)
        elif stage=='bpp':ret=serial.run_interface_bpp(fa,tmp,tag,cond)
        elif stage=='intarna':ret=serial.run_intarna(args,fa,tmp,tag,cond)
        else:raise ValueError(stage)
    finally:
        os.chdir(cwd)
    if Path(ret)!=table: raise RuntimeError('Unexpected serial output path')
    return read_table(table,stage,[p for p,_ in items])

def run_stage(stage,args,clean_fa,outdir,tag,condition,cofold_tsv=None):
    outdir=Path(outdir).resolve();clean_fa=Path(clean_fa).resolve();fa,table,manifest=paths(outdir,tag,condition,stage)
    if clean_fa.resolve()!=fa.resolve(): raise ValueError(f'Unexpected clean FASTA path: {clean_fa}')
    pairs=pairs_in(fa);fp=fingerprint(stage,args)
    cfs=read_table(cofold_tsv,'cofold',pairs) if stage=='ddg' else None
    existing=cache_load(outdir,tag,condition,stage,fp)
    if existing and cfs is not None and any(any(float(r[k])!=float(cfs[p][k]) for k in COLUMNS['cofold'][1:]) for p,r in existing[1].items()): existing=None
    if existing:
        print(f'[STRUCT] verified existing {stage}: {len(pairs)} rows',flush=True);return table
    # A result without a matching signature is not trusted.
    rows={};counts={}
    for d,c in sources(outdir,tag,condition,stage,fp):
        cache=cache_load(d,tag,c,stage,fp)
        if cache is None:continue
        ps,rs=cache;n=0
        for pid,row in rs.items():
            if pid in pairs and pid not in rows and pairs[pid]==ps[pid] and finite_row(row):
                rows[pid]=row;n+=1
        if n:counts[str(d)]=n
    todo=[(p,ab) for p,ab in pairs.items() if p not in rows]
    cpus=max(1,int(os.environ.get('SLURM_CPUS_PER_TASK','1')))
    workers=max(1,min(cpus,int(os.environ.get('DSRNA_STRUCT_WORKERS',str(cpus)))))
    chunk=max(1,int(os.environ.get('DSRNA_STRUCT_CHUNK','32')))
    # Ensure cached ddG's cofold inputs agree with the current cofold output.
    if cfs is not None:
        bad=[p for p,r in rows.items() if any(float(r[k])!=float(cfs[p][k]) for k in COLUMNS['cofold'][1:])]
        for p in bad:del rows[p];todo.append((p,pairs[p]))
    reused=len(rows);start=time.monotonic()
    print(f'[STRUCT] {stage} {condition}: total={len(pairs)} reused={reused} compute={len(todo)} workers={workers}',flush=True)
    if todo:
        # Each task invokes the unchanged historical routine on its own directory.
        # No shared output files or external-tool working directories across workers.
        with tempfile.TemporaryDirectory(prefix=f'dsrna_v2_{stage}_',dir=os.environ.get('TMPDIR')) as root:
            tasks=[]
            # Long sequences first distribute expensive work more evenly among workers.
            todo.sort(key=lambda x:len(x[1]['A'])+len(x[1]['B']),reverse=True)
            # Round-robin long/short candidates across chunks.
            nchunks=(len(todo)+chunk-1)//chunk
            for i in range(nchunks):
                items=todo[i::nchunks];cfrows=[cfs[p] for p,_ in items] if cfs is not None else []
                tasks.append((stage,args,items,cfrows,str(Path(root)/str(i)),tag,condition))
            done=0
            if workers==1:
                iterator=(work(t) for t in tasks)
                for r in iterator:
                    rows.update(r);done+=len(r)
                    if done==len(todo) or done% (chunk*4)==0:print(f'[STRUCT] {stage}: {done}/{len(todo)} computed; {time.monotonic()-start:.1f}s',flush=True)
            else:
                with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn')) as pool:
                    fs=[pool.submit(work,t) for t in tasks]
                    last=0
                    for future in as_completed(fs):
                        r=future.result();rows.update(r);done+=len(r)
                        if done==len(todo) or done-last>=chunk*4:
                            print(f'[STRUCT] {stage}: {done}/{len(todo)} computed; {time.monotonic()-start:.1f}s',flush=True);last=done
    if set(rows)!=set(pairs):raise RuntimeError('Result ID mismatch')
    write_table(table,stage,[rows[p] for p in pairs]);read_table(table,stage,pairs)
    mark_cache(outdir,tag,condition,stage,fp,{'reused':reused,'computed':len(todo),'sources':counts,'elapsed_seconds':time.monotonic()-start,'nonfinite_rows':sum(not finite_row(x) for x in rows.values())})
    return table

def run_rnacofold(args,clean_fa,outdir,tag,condition):
    return serial._load_pairs_from_clean_fasta(clean_fa),run_stage('cofold',args,clean_fa,outdir,tag,condition)

def run_ddg(args,clean_fa,cofold_tsv,outdir,tag,condition):
    return run_stage('ddg',args,clean_fa,outdir,tag,condition,cofold_tsv)

def run_interface_bpp(clean_fa,outdir,tag,condition):
    return run_stage('bpp',None,clean_fa,outdir,tag,condition)

def run_intarna(args,clean_fa,outdir,tag,condition):
    return run_stage('intarna',args,clean_fa,outdir,tag,condition)
