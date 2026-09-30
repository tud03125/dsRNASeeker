# dsRNASeeker Web Explorer tutorial

Public app: https://dsrnaseeker-kvqciapfzxrodulzmpdyb9.streamlit.app/

## ADPS example

1. Open **ADPS scoring**.
2. Select **Bundled example**.
3. Keep case=`CASE`, control=`CONTROL`, annotation policy=`conservative`.
4. Click **Calculate ADPS**.
5. Confirm a ranked table appears and download `dsRNASeeker_ranked.tsv`.

## Supervised example

1. Open **Same-study supervised**.
2. Select **Bundled example**.
3. Run grouped re-ranking.
4. Confirm fold diagnostics and supervised scores appear.

## Scope

The web interface is a lightweight candidate-table explorer.  Full FASTQ/BAM analysis,
alignment, TE quantification, RNA editing, rMATS, and candidate generation are provided
by the command-line dsRNASeeker workflow.

ADPS and supervised outputs are ranking scores, not calibrated probabilities of physical
duplex formation.
