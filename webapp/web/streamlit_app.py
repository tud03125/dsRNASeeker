from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from app.scoring import calculate_adps
from app.supervised import same_study_rerank
from app.registry import load_json

REPO_ROOT = ROOT.parent
EXAMPLES = ROOT / "examples"

st.set_page_config(
    page_title="dsRNASeeker Web Explorer",
    layout="wide",
)

st.title("dsRNASeeker Web Explorer")
st.caption(
    "Interactive candidate-table scoring and assay-matched supervised re-ranking "
    "for dsRNASeeker. The full FASTQ/BAM workflow is available through the command-line pipeline."
)

with st.sidebar:
    st.subheader("About")
    st.markdown(
        "dsRNASeeker prioritizes putative inverted transposable-element (TE) pairs "
        "using label-independent ADPS and optional assay-matched supervised re-ranking."
    )
    st.markdown(
        "[Source code on GitHub](https://github.com/tud03125/dsRNASeeker)"
    )
    st.caption(
        "ADPS is a ranking score, not a calibrated probability of physical dsRNA formation."
    )

t1, t2, t3, t4 = st.tabs(
    ["ADPS scoring", "Same-study supervised", "Benchmark registry", "Help & examples"]
)

def read_upload(f):
    return pd.read_csv(
        f,
        sep="\t" if f.name.endswith((".tsv", ".txt")) else ",",
        low_memory=False,
    )

with t1:
    st.subheader("Label-independent ADPS scoring")
    st.write(
        "Score a dsRNASeeker candidate/evidence table with the canonical adaptive ranking. "
        "Use the bundled example to verify formatting before uploading your own table."
    )
    source = st.radio(
        "Input source",
        ["Bundled example", "Upload file"],
        horizontal=True,
        key="adps_source",
    )
    df = None
    if source == "Bundled example":
        example = EXAMPLES / "example_adps_input.tsv"
        df = pd.read_csv(example, sep="\t")
        st.success(f"Loaded bundled example ({len(df)} candidate pairs).")
        st.download_button(
            "Download example ADPS input",
            example.read_bytes(),
            file_name="example_adps_input.tsv",
            mime="text/tab-separated-values",
        )
    else:
        f = st.file_uploader(
            "Upload dsRNASeeker candidate/evidence CSV or TSV",
            key="adps",
        )
        if f is not None:
            df = read_upload(f)

    case = st.text_input("Case label", value="CASE")
    control = st.text_input("Control label", value="CONTROL")
    policy = st.selectbox(
        "Annotation policy",
        ["conservative", "zrna_permissive"],
    )

    if df is not None:
        with st.expander("Preview input", expanded=False):
            st.dataframe(df.head(50), use_container_width=True)
        if st.button("Calculate ADPS", type="primary"):
            try:
                out = calculate_adps(df, case, control, policy)
                st.success(f"Scored {len(out)} candidate pairs.")
                cols = [
                    c for c in [
                        "pair_id", "adaptive_priority_score", "rank_score",
                        "priority_tier", "orientation_adps", "annotation_adps",
                        "case_expression_adps", "energy_adps", "interface_adps",
                        "case_editing_adps", "RI_adps",
                    ] if c in out.columns
                ]
                st.dataframe(out[cols].head(200), use_container_width=True)
                st.download_button(
                    "Download ranked TSV",
                    out.to_csv(sep="\t", index=False).encode(),
                    "dsRNASeeker_ranked.tsv",
                    mime="text/tab-separated-values",
                )
            except Exception as e:
                st.error(str(e))

with t2:
    st.subheader("Assay-matched supervised re-ranking")
    st.info(
        "Use candidate-level labels from the same assay/context. "
        "Current cross-study benchmarks do not support a universal pretrained classifier."
    )
    source = st.radio(
        "Input source",
        ["Bundled example", "Upload file"],
        horizontal=True,
        key="sup_source",
    )
    df = None
    if source == "Bundled example":
        example = EXAMPLES / "example_supervised_audit.tsv"
        df = pd.read_csv(example, sep="\t")
        st.success(f"Loaded bundled example ({len(df)} labeled candidate rows).")
        st.download_button(
            "Download example supervised input",
            example.read_bytes(),
            file_name="example_supervised_audit.tsv",
            mime="text/tab-separated-values",
        )
    else:
        f = st.file_uploader("Upload labeled audit TSV/CSV", key="sup")
        if f is not None:
            df = read_upload(f)

    folds = st.selectbox("Grouped outer folds", [3, 5], index=1)
    if df is not None:
        with st.expander("Preview input", expanded=False):
            st.dataframe(df.head(50), use_container_width=True)
        if st.button("Run grouped same-study re-ranking", type="primary"):
            try:
                ranked, fold_table = same_study_rerank(df, outer_folds=int(folds))
                st.success(f"Ranked {len(ranked)} labeled candidates.")
                st.markdown("**Fold assignment / diagnostic table**")
                if source == "Bundled example":
                    st.caption(
                        "The bundled 30-row example is a smoke-test fixture for the web workflow, "
                        "not a manuscript benchmark. Because individual grouped test folds are small, "
                        "fold-level AP/ROC-AUC values can be extreme (including 0 or 1) and should be "
                        "interpreted only as diagnostics of the example run."
                    )
                st.dataframe(fold_table, use_container_width=True)

                st.markdown("**Top supervised-ranked candidates**")
                preferred = [
                    "pair_id", "supervised_score", "label", "group_id",
                    "adaptive_priority_score", "orientation_adps", "annotation_adps",
                    "case_expression_adps", "energy_adps", "interface_adps",
                    "case_editing_adps", "RI_adps",
                ]
                display_cols = [c for c in preferred if c in ranked.columns]
                display_cols += [c for c in ranked.columns if c not in display_cols]
                st.dataframe(
                    ranked.loc[:, display_cols].head(200),
                    use_container_width=True,
                )
                st.download_button(
                    "Download supervised ranking TSV",
                    ranked.to_csv(sep="\t", index=False).encode(),
                    "dsRNASeeker_supervised_ranked.tsv",
                    mime="text/tab-separated-values",
                )
            except Exception as e:
                st.error(str(e))

with t3:
    st.subheader("Frozen benchmark context")
    studies = load_json("studies.json")
    models = load_json("models.json")

    st.markdown("**Benchmark studies**")
    study_df = pd.DataFrame(studies).copy()
    for c in ["primary_radius", "supervised_diagnostic_radius"]:
        if c not in study_df.columns:
            study_df[c] = None
    study_df = study_df[
        [
            "study", "assay", "primary_radius", "supervised_diagnostic_radius",
            "labels", "adps_ap", "adps_auc",
        ]
    ].rename(
        columns={
            "study": "Study",
            "assay": "Assay",
            "primary_radius": "Primary radius",
            "supervised_diagnostic_radius": "Supervised diagnostic radius",
            "labels": "Labels",
            "adps_ap": "ADPS AP",
            "adps_auc": "ADPS ROC-AUC",
        }
    ).fillna("—")
    st.dataframe(study_df, use_container_width=True, hide_index=True)

    st.caption(
        "A dash indicates that the radius is not designated for that role. "
        "For GSE184962/4, 5 kb is the expanded supervised-diagnostic radius, "
        "not the default primary radius."
    )

    st.markdown("**Ranking modes**")
    model_df = pd.DataFrame(models).rename(
        columns={
            "name": "Name",
            "type": "Type",
            "recommended_default": "Recommended default",
            "note": "Note",
        }
    )
    st.dataframe(model_df, use_container_width=True, hide_index=True)

    with st.expander("Raw registry JSON"):
        st.json(studies)
        st.json(models)

with t4:
    st.subheader("How to use this web interface")
    st.markdown(
        """
1. **ADPS scoring:** choose the bundled example, keep `CASE` and `CONTROL`,
   click **Calculate ADPS**, and inspect/download the ranked table.
2. **Same-study supervised:** choose the bundled labeled example and click
   **Run grouped same-study re-ranking**.
3. Replace the examples with your own dsRNASeeker candidate/evidence tables
   once the expected column structure is clear.
4. For full FASTQ/BAM processing, reference preparation, RNA-editing analysis,
   rMATS integration, and candidate generation, use the command-line workflow
   documented in the GitHub repository.

### Interpretation
- **ADPS** is label independent and intended for portable candidate ranking.
- **Supervised re-ranking** is assay matched; it should not be interpreted as a
  universal pretrained classifier.
- Structural and sequence-derived support are computational evidence, not direct
  proof of an in-vivo RNA duplex.
- The web explorer operates on candidate/evidence tables; it does not upload or
  expose private FCCC filesystem paths.

### Reproducibility
The GitHub repository provides installation instructions, environment specifications,
command-line examples, test inputs, and manuscript-reproduction materials. The exact
submission release will be identified in the manuscript by a versioned Git tag and a
matching permanent archive DOI once the release is frozen.
"""
    )
