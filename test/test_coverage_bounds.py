"""Dataset-independent smoke test for arm-aware chromosome-boundary clipping.
Run from the repository root using its tested production Conda env:
    python3 tests/test_coverage_bounds.py
It tests interval arithmetic on tiny synthetic bigWigs; no GEO data required.
"""
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import pandas as pd
import pyBigWig
from modules import coverage

def bwcreate(path, chrom_sizes):
    b=pyBigWig.open(str(path),'w')
    b.addHeader(list(chrom_sizes.items()))
    b.close()

def check():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        names=['s1.fwd.bw','s1.rev.bw']
        for n in names:bwcreate(root/n, {'chr1':1000,'chrAlt_example':200})
        samples=root/'samples.tsv';samples.write_text('sample_id\tcondition\ns1\tcase\n')
        args=SimpleNamespace(arm_pad=100,arm_min_cov=0.,samplesheet=str(samples),multibigwigsummary_exe='unused')
        pairs=pd.DataFrame([{'pair_id':'p1','A_chrom':'chr1','A_start':950,'A_end':995,'B_chrom':'chrAlt_example','B_start':5,'B_end':40}])
        records=[]
        def fake_multi(exe,bedfile,bw,out_npz,out_tsv):
            d=pd.read_csv(bedfile,sep='\t',header=None,names=['chr','start','end','name'])
            assert len(d)==2
            assert list(d['start'])==[850,0],d
            assert list(d['end'])==[1000,140],d
            records.append(tuple(d['name']))
            # Compatible test stub for deepTools BED-file raw-count format
            s=d[['chr','start','end']].copy();s['sample']=1.0
            s.to_csv(out_tsv,sep='\t',index=False)
        with patch.object(coverage,'_run_multi',side_effect=fake_multi):
            path=coverage.arm_aware_summaries(args,pairs,{'p1'},root,'inverted','case')
        assert len(records)==2
        got=pd.read_csv(path,sep='\t')
        assert len(got)==1 and got['pair_id'].iloc[0]=='p1'
        # Original unpadded arm cannot exceed chromosome end
        bad=pairs.copy();bad.loc[0,'A_end']=1001
        with patch.object(coverage,'_run_multi',side_effect=fake_multi):
            try:coverage.arm_aware_summaries(args,bad,{'p1'},root,'inverted','case')
            except ValueError as ex:assert 'Invalid unpadded TE arm' in str(ex)
            else:raise AssertionError('Bad unpadded arm was silently clipped')
        # A chromosome absent from the bigWig must not silently drop a pair
        bad=pairs.copy();bad.loc[0,'A_chrom']='chrMissing'
        with patch.object(coverage,'_run_multi',side_effect=fake_multi):
            try:coverage.arm_aware_summaries(args,bad,{'p1'},root,'inverted','case')
            except ValueError as ex:assert 'absent' in str(ex)
            else:raise AssertionError('Missing chromosome was silently ignored')
    print('PASS: terminal padding clipped; pair preserved; invalid originals rejected')
if __name__=='__main__':check()
