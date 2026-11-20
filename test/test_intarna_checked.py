import sys, tempfile, unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from modules import intarna_checked as m

class Checked(unittest.TestCase):
    def test_parse(self):
        for stdout,ret,expected in [('E\n-59.2\n',0,(-59.2,'SUCCESS')),('E\n',0,(None,'NO_INTERACTION'))]:
            with patch.object(m.subprocess,'run',return_value=SimpleNamespace(stdout=stdout,stderr='',returncode=ret)):
                self.assertEqual(m.predict('exe','AC','GU')[:2],expected)
        for stdout,ret in [('E\nNaN\n',0),('oops\n',0),('',1),('',0),('E\n1;2\n',0)]:
            with patch.object(m.subprocess,'run',return_value=SimpleNamespace(stdout=stdout,stderr='diagnostic',returncode=ret)):
                with self.assertRaises(RuntimeError):m.predict('exe','AC','GU')
    def test_cache_and_failure(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d); fa=d/'arms.fa';fa.write_text('>p|A\nAAAA\n>p|B\nUUUU\n>q|A\nCCCC\n>q|B\nGGGG\n')
            a=SimpleNamespace(intarna_exe='exe',threads=1); fingerprint={'exe':'exe','version':'test'}
            out=d/'duplex_pairs.inv.case.IntaRNA.tsv';out.write_text('pair_id\tE\np\t\nq\t\n')
            with patch.object(m,'preflight_intarna',return_value=fingerprint),patch.object(m,'predict',side_effect=[(-1.,'SUCCESS',''),(None,'NO_INTERACTION','')]) as f:
                self.assertEqual(m.run_intarna_checked(a,fa,d,'inv','case'),out);self.assertEqual(f.call_count,2)
                self.assertEqual(m.run_intarna_checked(a,fa,d,'inv','case'),out);self.assertEqual(f.call_count,2)
            fa.write_text(fa.read_text().replace('AAAA','AAAC'))
            with patch.object(m,'preflight_intarna',return_value=fingerprint),patch.object(m,'predict',side_effect=RuntimeError('bad')):
                with self.assertRaises(RuntimeError):m.run_intarna_checked(a,fa,d,'inv','case')
            self.assertTrue(out.with_suffix('.failure.json').exists()); self.assertFalse(m.reuse_ok(out,out.with_suffix('.status.tsv'),out.with_suffix('.provenance.json'),fingerprint,fa,['p','q']))
    def test_invalid_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';p.write_text('>p|A\nAAA\n')
            with self.assertRaises(ValueError):m.load_pairs(p)
if __name__=='__main__':unittest.main()
