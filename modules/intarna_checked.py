"""Checked, ancillary IntaRNA predictions. No change to ADPS or folding models."""
from __future__ import annotations
import csv, hashlib, io, json, math, os, shutil, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def atomic(path, content):
    path=Path(path)
    with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,delete=False,encoding='utf-8') as f:
        f.write(content); tmp=Path(f.name)
    tmp.replace(path)


def predict(exe, a, b, timeout=None):
    cmd=[exe,'--query',a,'--target',b,'--outMode','C','--outCsvCols','E']
    p=subprocess.run(cmd,capture_output=True,text=True,check=False,timeout=timeout)
    if p.returncode:
        raise RuntimeError(f'IntaRNA EXECUTION_ERROR exit={p.returncode}: {p.stderr[-2000:]}')
    rows=list(csv.reader(io.StringIO(p.stdout),delimiter=';'))
    rows=[r for r in rows if any(x.strip() for x in r)]
    if not rows or rows[0]!=['E']:
        raise RuntimeError(f'IntaRNA PARSE_ERROR: expected E header; stdout={p.stdout[:500]!r}; stderr={p.stderr[:500]!r}')
    if len(rows)==1:
        return None,'NO_INTERACTION',p.stderr[-2000:]
    values=[]
    for r in rows[1:]:
        try:
            if len(r)!=1: raise ValueError(r)
            v=float(r[0]); assert math.isfinite(v)
        except (ValueError,AssertionError) as exc:
            raise RuntimeError(f'IntaRNA PARSE_ERROR: {r!r}') from exc
        values.append(v)
    # Preserve the historical last-line energy convention, not a new scoring formula.
    return values[-1],'SUCCESS',p.stderr[-2000:]


@lru_cache(maxsize=8)
def _probe(exe, digest):
    p=subprocess.run([exe,'--version'],capture_output=True,text=True,timeout=120)
    if p.returncode: raise RuntimeError(f'IntaRNA preflight: incompatible executable {exe}: {p.stderr}')
    q='GGGAGGAGGGAGGAGGGAGGAGGG'; t=q.translate(str.maketrans('ACGU','UGCA'))[::-1]
    e,status,_=predict(exe,q,t,timeout=120)
    if status!='SUCCESS': raise RuntimeError(f'IntaRNA preflight returned no finite energy: {exe}')
    return {'exe':exe,'exe_sha256':digest,'version':p.stdout.strip() or p.stderr.strip(), 'probe_energy':e,
            'module_sha256':sha(__file__), 'python':sys.version, 'ld_library_path':os.environ.get('LD_LIBRARY_PATH','')}


def preflight_intarna(exe):
    p=shutil.which(exe)
    if not p: raise RuntimeError(f'IntaRNA not found: {exe}; use --intarna-exe with a tested modern executable')
    p=str(Path(p).resolve())
    return _probe(p,sha(p))


def load_pairs(fa):
    pairs={}; seen=set(); cur=None
    for ln in Path(fa).read_text().splitlines():
        if ln.startswith('>'):
            pid,arm=ln[1:].strip().rsplit('|',1)
            if not pid or arm not in ('A','B') or (pid,arm) in seen: raise ValueError(f'Invalid duplex FASTA header {ln}')
            cur=(pid,arm); seen.add(cur); pairs.setdefault(pid,{})[arm]=''
        elif ln.strip():
            if cur is None: raise ValueError('Sequence before header')
            pairs[cur[0]][cur[1]]+=ln.strip()
    if any(set(x)!={'A','B'} or not x['A'] or not x['B'] for x in pairs.values()): raise ValueError('Incomplete duplex arms')
    return pairs


def tsv(header, rows):
    b=io.StringIO(); w=csv.writer(b,delimiter='\t',lineterminator='\n'); w.writerow(header); w.writerows(rows); return b.getvalue()


def reuse_ok(out, status, meta, fp, fa, ids):
    try:
        d=json.loads(meta.read_text())
        if d['fingerprint']!=fp or d['fasta_sha256']!=sha(fa) or d['table_sha256']!=sha(out) or d['status_sha256']!=sha(status): return False
        with out.open() as f: r=list(csv.DictReader(f,delimiter='\t'))
        with status.open() as f: s=list(csv.DictReader(f,delimiter='\t'))
        if [x['pair_id'] for x in r]!=ids or [x['pair_id'] for x in s]!=ids: return False
        return all((y['status']=='SUCCESS' and math.isfinite(float(x['E']))) or (y['status']=='NO_INTERACTION' and x['E']=='') for x,y in zip(r,s))
    except (OSError,ValueError,KeyError,TypeError): return False


def run_intarna_checked(args, clean_fa, outdir, tag, condition):
    fp=preflight_intarna(args.intarna_exe); pairs=load_pairs(clean_fa); ids=list(pairs)
    out=Path(outdir)/f'duplex_pairs.{tag}.{condition}.IntaRNA.tsv'
    status=out.with_suffix('.status.tsv'); meta=out.with_suffix('.provenance.json')
    if reuse_ok(out,status,meta,fp,clean_fa,ids):
        print(f'[IntaRNA] verified cache: {len(ids)} pairs',flush=True);return out
    # Invalidate prior manifest before recomputation. Preserve old TSV for diagnosis only.
    meta.unlink(missing_ok=True)
    # Old tables are never accepted without validated provenance. Do not publish a new partial table.
    workers=max(1,min(int(getattr(args,'threads',1)),int(os.environ.get('SLURM_CPUS_PER_TASK',getattr(args,'threads',1))),int(os.environ.get('DSRNA_INTARNA_WORKERS','1'))))
    def one(pid):
        ab=pairs[pid]
        try: return predict(fp['exe'],ab['A'],ab['B'])
        except Exception as e: raise RuntimeError(f'IntaRNA failed for {pid}: {e}') from e
    results=[]
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            # Bounded submission: a failure cannot leave thousands of pending predictions.
            for offset in range(0,len(ids),max(1,workers*2)):
                for x in pool.map(one,ids[offset:offset+max(1,workers*2)]):
                    results.append(x)
                    i=len(results)
                    if i%100==0 or i==len(ids): print(f'[IntaRNA] {i}/{len(ids)}',flush=True)
    except Exception as exc:
        atomic(out.with_suffix('.failure.json'),json.dumps({'error':str(exc),'fingerprint':fp,'fasta_sha256':sha(clean_fa)},indent=2)+'\n')
        raise
    atomic(out,tsv(['pair_id','E'],[(pid,'' if x[0] is None else x[0]) for pid,x in zip(ids,results)]))
    atomic(status,tsv(['pair_id','status','stderr'],[(pid,x[1],x[2]) for pid,x in zip(ids,results)]))
    atomic(meta,json.dumps({'fingerprint':fp,'fasta_sha256':sha(clean_fa),'table_sha256':sha(out),'status_sha256':sha(status), 'n_success':sum(x[1]=='SUCCESS' for x in results),'n_no_interaction':sum(x[1]=='NO_INTERACTION' for x in results)},indent=2)+'\n')
    # Remove stale failure report only after successful complete output.
    out.with_suffix('.failure.json').unlink(missing_ok=True)
    return out
