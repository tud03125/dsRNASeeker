from __future__ import annotations
from pathlib import Path
import re, random, subprocess, hashlib
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
from .utils import run_cmd


def _expected_complete_pair_count(clean_fa):
    """
    Count duplex pairs that have both non-empty A and B arms in the clean FASTA.
    This is used to decide whether an existing output TSV is complete enough to skip.
    """
    pairs = _load_pairs_from_clean_fasta(clean_fa)
    return sum(
        1
        for ab in pairs.values()
        if "A" in ab and "B" in ab and bool(ab["A"]) and bool(ab["B"])
    )


def _valid_tsv(path, required_cols=None, min_rows=1):
    """
    Return True only if a TSV exists, is readable, has the expected columns,
    and has at least min_rows data rows.

    Important: this prevents a failed/partial/empty output from being treated
    as complete. The write functions below write to *.tmp first and then rename
    atomically, which further reduces the risk of partial files on reruns.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return False

    try:
        # For large files, do not load the full table if we only need to verify
        # existence/columns/row-count threshold. nrows=0 is allowed for min_rows=0.
        if min_rows and int(min_rows) > 0:
            df = pd.read_csv(path, sep="\t", nrows=int(min_rows))
        else:
            df = pd.read_csv(path, sep="\t", nrows=0)
    except Exception:
        return False

    if required_cols:
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            return False

    if min_rows and int(min_rows) > 0:
        return len(df) >= int(min_rows)

    return True


def _write_tsv_atomic(df, out):
    """Write TSV through a temporary file, then atomically replace the target."""
    out = Path(out)
    tmp = out.with_suffix(out.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    tmp.replace(out)



def prepare_duplex_inputs(bedtools_exe, fasta, pairs, kept_pair_ids, outdir, tag, condition):
    """
    Extract the two TE arms in a common genomic orientation.

    A_strand/B_strand describe RepeatMasker insertion orientation. They are
    appropriate for deciding whether a pair is inverted, but must not be used
    to reverse-complement the two arms independently for duplex thermodynamics.
    """
    outdir = Path(outdir)
    pairs = pairs[pairs['pair_id'].isin(kept_pair_ids)]
    rows = []
    for _, r in pairs.iterrows():
        rows.append([r['A_chrom'], r['A_start'], r['A_end'], f"{r['pair_id']}|A"])
        rows.append([r['B_chrom'], r['B_start'], r['B_end'], f"{r['pair_id']}|B"])

    bed = outdir / f'duplex_arms.{tag}.{condition}.bed'
    pd.DataFrame(rows, columns=['chrom','start','end','name']).to_csv(
        bed, sep='\t', header=False, index=False
    )
    cp = run_cmd(
        [bedtools_exe, 'getfasta', '-fi', str(fasta), '-bed', str(bed), '-name'],
        capture=True,
        cwd=str(outdir),
    )
    fa = outdir / f'duplex_arms.{tag}.{condition}.fa'
    fa.write_text(cp.stdout)

    clean = outdir / f'duplex_arms.{tag}.{condition}.clean.fa'
    with clean.open('w') as fo, fa.open() as fi:
        for line in fi:
            if line.startswith('>'):
                h = line[1:].strip().split()[0].split('::', 1)[0]
                if '|' not in h:
                    continue
                pair, arm = h.split('|', 1)
                arm = arm[:1]
                if arm not in ('A', 'B'):
                    continue
                fo.write(f'>{pair}|{arm}\n')
            else:
                fo.write(line.upper().replace('T', 'U'))

    provenance = outdir / f'duplex_arms.{tag}.{condition}.sequence_orientation.txt'
    provenance.write_text(
        "sequence_orientation=common_reference_plus\n"
        "repeatmasker_strands_used_for_pair_orientation_not_independent_sequence_reverse_complement\n"
    )
    return clean

def _load_pairs_from_clean_fasta(clean_fa):
    fa = Path(clean_fa).read_text().splitlines()
    pairs={}; pid=arm=None
    for line in fa:
        if line.startswith('>'):
            pid, arm = line[1:].strip().split('|',1)
            pairs.setdefault(pid,{})[arm] = ''
        else:
            if pid and arm:
                pairs[pid][arm] += line.strip()
    return pairs


def run_rnacofold(args, clean_fa, outdir, tag, condition):
    outdir=Path(outdir)
    pairs=_load_pairs_from_clean_fasta(clean_fa)
    expected_n = _expected_complete_pair_count(clean_fa)

    out_tsv=outdir/f'duplex_pairs.{tag}.{condition}.cofold_mfe.tsv'
    required_cols = ['pair_id','RNAcofold_MFE_kcalmol','lenA','lenB','len_total','MFE_norm_kcalpermkb']
    if _valid_tsv(out_tsv, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid RNAcofold MFE file found: {out_tsv}")
        return pairs, out_tsv

    cofold_in=outdir/f'duplex_pairs.{tag}.{condition}.cofold.in'
    lengths=outdir/f'duplex_pairs.{tag}.{condition}.lengths.tsv'
    with cofold_in.open('w') as fo, lengths.open('w') as fl:
        fl.write('pair_id\tlenA\tlenB\tlen_total\n')
        for pid,ab in pairs.items():
            if 'A' in ab and 'B' in ab and ab['A'] and ab['B']:
                fo.write(f'>{pid}\n{ab["A"]}&{ab["B"]}\n')
                fl.write(f'{pid}\t{len(ab["A"])}\t{len(ab["B"])}\t{len(ab["A"])+len(ab["B"])}\n')
    with cofold_in.open("r") as fin:
      proc = subprocess.run(
          [args.rnacofold_exe, "--noPS"],
          stdin=fin,
          text=True,
          capture_output=True,
          check=True,
          cwd=outdir,
      )
    (outdir/f'duplex_pairs.{tag}.{condition}.cofold').write_text(proc.stdout)
    lines=proc.stdout.splitlines(); out=[]; cur=None
    for i,line in enumerate(lines):
        if line.startswith('>'): cur=line[1:].strip()
        elif cur:
            m = re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', line) or (re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', lines[i+1]) if i+1 < len(lines) else None)
            if m:
                out.append((cur, float(m.group(1)))); cur=None
    mfe=pd.DataFrame(out, columns=['pair_id','RNAcofold_MFE_kcalmol'])
    lens=pd.read_csv(lengths, sep='\t')
    X=mfe.merge(lens,on='pair_id', how='left')
    X['MFE_norm_kcalpermkb']=X['RNAcofold_MFE_kcalmol']/(X['len_total']/1000.0)
    _write_tsv_atomic(X, out_tsv)
    return pairs, out_tsv


def run_ddg(args, clean_fa, cofold_tsv, outdir, tag, condition):
    outdir=Path(outdir)
    expected_n = _expected_complete_pair_count(clean_fa)
    out=outdir/f'duplex_pairs.{tag}.{condition}.ddg.tsv'
    required_cols = [
        'pair_id','RNAcofold_MFE_kcalmol','lenA','lenB','len_total',
        'MFE_norm_kcalpermkb','RNAfold_A_MFE_kcalmol','RNAfold_B_MFE_kcalmol',
        'ddG_interaction_kcalmol','ddG_norm_kcalpermkb'
    ]
    if _valid_tsv(out, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid ddG file found: {out}")
        return out

    pairs=_load_pairs_from_clean_fasta(clean_fa)
    Ain=outdir/f'duplex_pairs.{tag}.{condition}.A.fold.in'; Bin=outdir/f'duplex_pairs.{tag}.{condition}.B.fold.in'
    with Ain.open('w') as fA, Bin.open('w') as fB:
        for pid,ab in pairs.items():
            if 'A' in ab and 'B' in ab and ab['A'] and ab['B']:
                fA.write(f'>{pid}\n{ab["A"]}\n'); fB.write(f'>{pid}\n{ab["B"]}\n')
    with Ain.open("r") as finA:
        outA = subprocess.run(
            [args.rnafold_exe, "--noPS"],
            stdin=finA,
            text=True,
            capture_output=True,
            check=True,
            cwd=outdir,
        ).stdout

    with Bin.open("r") as finB:
        outB = subprocess.run(
            [args.rnafold_exe, "--noPS"],
            stdin=finB,
            text=True,
            capture_output=True,
            check=True,
            cwd=outdir,
        ).stdout
    def parse_fold(text):
        lines=text.splitlines(); out=[]; cur=None
        for ln in lines:
            if ln.startswith('>'): cur=ln[1:].strip()
            else:
                m=re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', ln)
                if m and cur:
                    out.append((cur, float(m.group(1)))); cur=None
        return pd.DataFrame(out, columns=['pair_id','G_kcalmol'])
    A=parse_fold(outA).rename(columns={'G_kcalmol':'RNAfold_A_MFE_kcalmol'})
    B=parse_fold(outB).rename(columns={'G_kcalmol':'RNAfold_B_MFE_kcalmol'})
    C=pd.read_csv(cofold_tsv, sep='\t')
    X=C.merge(A,on='pair_id',how='left').merge(B,on='pair_id',how='left')
    X['ddG_interaction_kcalmol']=X['RNAcofold_MFE_kcalmol']-(X['RNAfold_A_MFE_kcalmol']+X['RNAfold_B_MFE_kcalmol'])
    X['ddG_norm_kcalpermkb']=X['ddG_interaction_kcalmol']/(X['len_total']/1000.0)
    _write_tsv_atomic(X, out)
    return out



def run_interface_bpp(clean_fa, outdir, tag, condition):
    """
    Summarize predicted cross-arm base-pair probabilities for the isolated A+B
    cofold model. This is not an in-vivo encounter probability.
    """
    import RNA
    outdir = Path(outdir)
    expected_n = _expected_complete_pair_count(clean_fa)
    out = outdir / f'duplex_pairs.{tag}.{condition}.interface_bpp.tsv'
    required_cols = [
        'pair_id', 'interface_bpp_sum', 'interface_bpp_max', 'interface_bpp_n',
        'interface_bpp_n_ge_1e5',
        'interface_bpp_expected_fraction_shorter',
        'interface_bpp_mean_arm_fraction',
    ]
    if _valid_tsv(out, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid ViennaRNA interface BPP file found: {out}")
        return out

    pairs = _load_pairs_from_clean_fasta(clean_fa)
    rows = []
    for pid, ab in pairs.items():
        if 'A' not in ab or 'B' not in ab:
            continue
        A, B = ab['A'], ab['B']
        if not A or not B:
            continue

        seq = A + '&' + B
        lenA, lenB = len(A), len(B)
        n = lenA + lenB
        fc = RNA.fold_compound(seq)
        if hasattr(fc, "pf_dimer"):
            fc.pf_dimer()
        else:
            fc.pf()
        bppm = fc.bpp()
        one_based = (len(bppm) == n + 1)

        def get_p(i, j):
            return float(bppm[i][j]) if one_based else float(bppm[i-1][j-1])

        p_sum = 0.0
        max_p = 0.0
        p_n = 0
        p_n_ge_1e5 = 0
        for i in range(1, lenA + 1):
            for j in range(lenA + 1, n + 1):
                p = get_p(i, j)
                if p > 0:
                    p_sum += p
                    p_n += 1
                    if p >= 1e-5:
                        p_n_ge_1e5 += 1
                    max_p = max(max_p, p)

        shorter = float(min(lenA, lenB))
        expected_fraction_shorter = p_sum / shorter if shorter > 0 else np.nan
        mean_arm_fraction = (
            0.5 * ((p_sum / float(lenA)) + (p_sum / float(lenB)))
            if lenA > 0 and lenB > 0 else np.nan
        )

        rows.append({
            'pair_id': pid,
            'interface_bpp_sum': float(p_sum),
            'interface_bpp_max': float(max_p),
            'interface_bpp_n': int(p_n),
            'interface_bpp_n_ge_1e5': int(p_n_ge_1e5),
            'interface_bpp_expected_fraction_shorter': float(expected_fraction_shorter),
            'interface_bpp_mean_arm_fraction': float(mean_arm_fraction),
        })

    _write_tsv_atomic(pd.DataFrame(rows), out)
    return out


def _dinuc_counts(seq: str) -> Counter:
    """Return exact adjacent dinucleotide counts for an RNA/DNA sequence."""
    return Counter(zip(seq[:-1], seq[1:]))


def _dinuc_shuffle_exact(seq: str, rng: random.Random) -> str:
    """
    Randomize a sequence while preserving its exact dinucleotide multiset.

    The sequence is represented as a directed multigraph of adjacent bases and
    an Eulerian trail is generated after randomizing outgoing edge order.  The
    resulting permutation is verified before it is returned.
    """
    if len(seq) < 2:
        return seq
    edges = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        edges[a].append(b)
    for key in list(edges):
        rng.shuffle(edges[key])

    stack = [seq[0]]
    path = []
    while stack:
        vertex = stack[-1]
        if edges[vertex]:
            stack.append(edges[vertex].pop())
        else:
            path.append(stack.pop())
    shuffled = ''.join(reversed(path))
    if len(shuffled) != len(seq) or _dinuc_counts(shuffled) != _dinuc_counts(seq):
        raise RuntimeError("Exact dinucleotide shuffle invariant failed")
    return shuffled


def run_null_z(args, clean_fa, ddg_tsv, outdir, tag, condition):
    """
    Exact dinucleotide-preserving null for ddG. Shuffled monomer MFEs are
    recomputed for every null sequence before ddG is calculated.
    """
    outdir = Path(outdir)
    expected_n = _expected_complete_pair_count(clean_fa)
    out = outdir / f'duplex_pairs.{tag}.{condition}.nullZ.tsv'
    required_cols = [
        'pair_id','ddG_mu_null','ddG_sd_null','ddG_Z',
        'null_n_requested','null_n_effective','null_shuffle_exact_dinuc','null_shuffle_method'
    ]
    if _valid_tsv(out, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid null-Z file found: {out}")
        return out

    pairs = _load_pairs_from_clean_fasta(clean_fa)
    obs = pd.read_csv(ddg_tsv, sep='\t')
    def run_batch(exe, input_path):
        with input_path.open("r") as fin:
            proc = subprocess.run(
                [exe, "--noPS"], stdin=fin, text=True, capture_output=True,
                check=True, cwd=outdir
            )
        vals, cur = {}, None
        lines = proc.stdout.splitlines()
        for i, ln in enumerate(lines):
            if ln.startswith('>'):
                cur = ln[1:].strip()
                continue
            if cur:
                m = re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', ln)
                if not m and i + 1 < len(lines):
                    m = re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', lines[i+1])
                if m:
                    vals[cur] = float(m.group(1))
                    cur = None
        return vals

    rows = []
    for _, r in obs.iterrows():
        pid = str(r['pair_id'])
        A = pairs.get(pid, {}).get('A', '')
        B = pairs.get(pid, {}).get('B', '')
        if not A or not B:
            continue

        nnull = int(args.null_n)
        # Stable per-pair RNG makes the null reproducible even if candidate row
        # order changes between runs.
        seed_material = f"{int(args.null_seed)}::{pid}".encode('utf-8')
        pair_seed = int.from_bytes(hashlib.sha256(seed_material).digest()[:8], 'big')
        rng = random.Random(pair_seed)
        shuffled = [(k, _dinuc_shuffle_exact(A, rng), _dinuc_shuffle_exact(B, rng)) for k in range(nnull)]
        token = hashlib.md5(pid.encode('utf-8')).hexdigest()[:12]
        co_in = outdir / f'.null_{token}.cofold.in'
        a_in = outdir / f'.null_{token}.A.fold.in'
        b_in = outdir / f'.null_{token}.B.fold.in'
        with co_in.open('w') as fc, a_in.open('w') as fa, b_in.open('w') as fb:
            for k, As, Bs in shuffled:
                name = f'{token}__null{k}'
                fc.write(f'>{name}\n{As}&{Bs}\n')
                fa.write(f'>{name}\n{As}\n')
                fb.write(f'>{name}\n{Bs}\n')

        co = run_batch(args.rnacofold_exe, co_in)
        am = run_batch(args.rnafold_exe, a_in)
        bm = run_batch(args.rnafold_exe, b_in)

        vals = []
        for k, _, _ in shuffled:
            name = f'{token}__null{k}'
            if name in co and name in am and name in bm:
                vals.append(co[name] - (am[name] + bm[name]))
        arr = np.asarray(vals, dtype=float)
        mu = float(np.nanmean(arr)) if len(arr) else np.nan
        sd = float(np.nanstd(arr, ddof=1)) if len(arr) > 1 else np.nan
        z = ((float(r['ddG_interaction_kcalmol']) - mu) / sd
             if pd.notna(sd) and sd > 0 else np.nan)
        rows.append({
            'pair_id': pid, 'ddG_mu_null': mu, 'ddG_sd_null': sd, 'ddG_Z': z,
            'null_n_requested': nnull, 'null_n_effective': int(len(arr)),
            'null_shuffle_exact_dinuc': True, 'null_shuffle_method': 'randomized_eulerian_exact_dinucleotide',
        })
        for p in (co_in, a_in, b_in):
            try:
                p.unlink()
            except FileNotFoundError:
                pass

    _write_tsv_atomic(pd.DataFrame(rows), out)
    return out

def run_intarna(args, clean_fa, outdir, tag, condition):
    outdir=Path(outdir)
    expected_n = _expected_complete_pair_count(clean_fa)
    out=outdir/f'duplex_pairs.{tag}.{condition}.IntaRNA.tsv'
    required_cols = ['pair_id','E']
    if _valid_tsv(out, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid IntaRNA file found: {out}")
        return out

    pairs=_load_pairs_from_clean_fasta(clean_fa); rows=[]
    for pid,ab in pairs.items():
        if 'A' not in ab or 'B' not in ab: continue
        try:
            proc=subprocess.run([args.intarna_exe, '--query', ab['A'], '--target', ab['B'], '--outMode', 'C', '--outCsvCols', 'E'], text=True, capture_output=True, check=True)
            val=proc.stdout.strip().splitlines()[-1].strip()
            energy=float(val) if val not in ('','E') else np.nan
        except Exception:
            energy=np.nan
        rows.append({'pair_id':pid,'E':energy})
    _write_tsv_atomic(pd.DataFrame(rows), out)
    return out
