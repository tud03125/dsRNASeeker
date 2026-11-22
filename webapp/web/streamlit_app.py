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

def read_upload(f):
    return pd.read_csv(
        f,
        sep="\t" if f.name.endswith((".tsv", ".txt")) else ",",
        low_memory=False,
    )


def adps_page():
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


def supervised_page():
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



# PORTAL_REVISION_20261118

BASE = (
    REPO_ROOT / "reproduce_manuscript"
    / "benchmark_sensitivity"
)

GITHUB = "https://github.com/tud03125/dsRNASeeker"

SOFTWARE_DOI = "https://doi.org/10.5281/zenodo.23407322"

DATA_DOI = "https://doi.org/10.5281/zenodo.24296933"

WEBSITE = "https://dsrnaseeker.streamlit.app/"

EXAMPLE_PAIR = (
    "IAPLTR1a_Mm_dup1094__RLTR10C_dup1140"
)

PAGES = [
    "Home",
    "Candidate explorer",
    "Run an analysis",
    "Benchmarks",
    "Documentation",
    "Citation",
]

CHANNELS = {
    "orientation_adps": "Pair orientation",
    "annotation_adps": "TE annotation",
    "case_expression_adps": "Case-associated TE expression",
    "energy_adps": "Duplex energetics",
    "interface_adps": "Pairing interface",
    "case_editing_adps": "A-to-I editing",
    "RI_adps": "Intron retention",
}

STUDIES = {
    "FLAG_308489": "FLAG-ZBP1",
    "J2_244103": "J2",
    "Z22_184962": "Z22 (diagnostic)",
}

@st.cache_data
def load_table(relative_path):
    return pd.read_csv(BASE / relative_path)


def navigate(page):
    st.session_state["portal_page"] = page


def home_page():
    st.header("Discover candidate endogenous dsRNAs")

    st.write(
        "dsRNASeeker integrates RNA-seq evidence "
        "to prioritize transposable-element pairs "
        "that may form endogenous double-stranded RNA."
    )

    c1, c2 = st.columns(2)

    with c1:
        st.button(
            "Explore an example",
            key="home_example",
            type="primary",
            on_click=navigate,
            args=("Candidate explorer",),
        )

    with c2:
        st.button(
            "Run an analysis",
            key="home_analysis",
            on_click=navigate,
            args=("Run an analysis",),
        )

    st.divider()

    st.subheader("How does dsRNASeeker work?")

    c1, c2, c3 = st.columns(3)

    with c1:
        st.markdown("**1. Candidate discovery**")
        st.write(
            "Identify nearby transposable-element pairs "
            "from RNA-seq and genome annotations."
        )

    with c2:
        st.markdown("**2. Evidence integration**")
        st.write(
            "Combine TE expression, genomic context, "
            "RNA pairing, editing and splicing evidence."
        )

    with c3:
        st.markdown("**3. Candidate ranking**")
        st.write(
            "Prioritize candidates using ADPS, with "
            "optional assay-matched supervised ranking."
        )

    st.info(
        "The web explorer analyzes processed candidate "
        "tables. The full FASTQ/BAM workflow is available "
        "through the command-line pipeline."
    )

    c1, c2 = st.columns(2)

    with c1:
        st.button(
            "Documentation",
            key="home_docs",
            on_click=navigate,
            args=("Documentation",),
        )

    with c2:
        st.button(
            "Cite dsRNASeeker",
            key="home_cite",
            on_click=navigate,
            args=("Citation",),
        )


def candidate_page():
    st.header("Candidate explorer")

    st.caption(
        "FLAG-ZBP1 benchmark | GSE308488/9 | "
        "HSV-1 versus mock | mm39 | 1 kb"
    )

    d = load_table(
        "01_benchmarks/"
        "FLAG_308489_verified_labeled_scores.csv"
    )

    mapping = load_table("03_transfer/ARM_MAPPING.csv")

    mapping = mapping.loc[
        mapping["study"].eq("FLAG_308489")
        & mapping["status"].eq("OK")
    ].copy()

    # Offer only verified inverted positives for
    # which both arm-coordinate records are available.
    complete = (
        mapping.groupby("pair_id")["arm"]
        .agg(lambda x: len(x) == 2 and set(x) == {"A", "B"})
    )

    eligible = d.loc[
        d.label.eq(1)
        & d.genomic_orientation.eq("inverted")
        & d.pair_id.isin(complete[complete].index)
    ].copy()

    eligible = eligible.sort_values(
        "ADPS", ascending=False
    )

    options = eligible.pair_id.tolist()

    assert EXAMPLE_PAIR in options

    pair = st.selectbox(
        "Select a manuscript benchmark candidate",
        options,
        index=options.index(EXAMPLE_PAIR),
        key="example_pair",
    )

    r = eligible.loc[eligible.pair_id.eq(pair)].iloc[0]

    a = mapping.loc[mapping.pair_id.eq(pair)].copy()

    assert len(a) == 2
    assert set(a["arm"]) == {"A", "B"}

    st.subheader(
        f"Gene-group annotation: {r['group_id']}"
    )

    c1, c2, c3 = st.columns(3)

    c1.metric("ADPS", f"{r['ADPS']:.3f}")
    c2.metric("Supervised score", f"{r['Supervised']:.3f}")
    c3.metric("Assay label", "Positive")

    st.caption(
        "The supervised score is an archived held-out "
        "prediction. The example is a candidate from "
        "the paper's benchmark, not a newly predicted pair."
    )

    st.subheader("Transposable-element arms")

    a["TE locus"] = a["arm"].map({
        "A": r["A_TE_id"],
        "B": r["B_TE_id"],
    })

    coordinates = a["coordinates"].str.extract(
        r"^([^:]+):(\d+):(\d+)$"
    )

    if coordinates.notna().all().all():
        a["Chromosome"] = coordinates[0]
        a["Start"] = coordinates[1].astype(int)
        a["End"] = coordinates[2].astype(int)

        a = a.sort_values(["Chromosome", "Start"])

    c1, c2 = st.columns(2)

    for col, (_, arm) in zip((c1, c2), a.iterrows()):
        with col:
            st.markdown(f"**Arm {arm['arm']}**")
            st.write(arm["TE locus"])
            st.code(arm["coordinates"], language=None)
            st.metric(
                "Reference-arm length",
                f"{int(arm['length'])} nt"
            )

    st.markdown("**Genomic orientation: inverted**")

    st.caption(
        "Arms are shown in genomic-coordinate order "
        "where parseable. Coordinates are reproduced "
        "as archived in the source mapping."
    )

    st.dataframe(
        a[["arm", "TE locus", "coordinates", "length"]],
        hide_index=True,
        use_container_width=True,
    )

    st.subheader("Supporting molecular evidence")

    evidence = pd.DataFrame({
        "Evidence channel": list(CHANNELS.values()),
        "Evidence score": [
            float(r[k]) for k in CHANNELS
        ],
    })

    st.bar_chart(
        evidence.set_index("Evidence channel")
    )

    st.dataframe(
        evidence.round(3),
        hide_index=True,
        use_container_width=True,
    )

    st.caption(
        "These are normalized evidence-component "
        "scores used for ranking. They are not raw "
        "expression fold changes, editing percentages "
        "or probabilities of physical duplex formation."
    )

    st.download_button(
        "Download this candidate record",
        r.to_frame().T.to_csv(index=False),
        file_name="dsRNASeeker_benchmark_candidate.csv",
        mime="text/csv",
    )

    st.markdown(
        "[View GEO study]"
        "(https://www.ncbi.nlm.nih.gov/geo/query/"
        "acc.cgi?acc=GSE308489)"
        " | "
        f"[Archived benchmark data]({DATA_DOI})"
    )


def analysis_page():
    st.header("Run an analysis")

    st.info(
        "The online analysis accepts processed candidate "
        "and evidence tables. For complete analysis from "
        "FASTQ or BAM, use the command-line workflow."
    )

    st.warning(
        "Do not upload confidential or identifiable "
        "patient data to this public cloud application."
    )

    mode = st.radio(
        "Analysis type",
        [
            "Label-independent ADPS scoring",
            "Assay-matched supervised ranking",
        ],
        key="analysis_type",
    )

    if mode == "Label-independent ADPS scoring":
        adps_page()
    else:
        supervised_page()

    st.markdown(
        f"[Full FASTQ/BAM workflow]({GITHUB})"
    )


def benchmark_page():
    st.header("Benchmark performance")

    st.write(
        "Archived performance from the manuscript's "
        "assay-labeled benchmark datasets."
    )

    d = load_table(
        "01_benchmarks/"
        "ALL_CHANNELS_AND_STRICT_SUBSET.csv"
    )

    subset = st.radio(
        "Candidate subset",
        ["All labeled candidates",
         "Strictly inverted labeled candidates"],
        horizontal=True,
        key="benchmark_subset",
    )

    key = (
        "all" if subset == "All labeled candidates"
        else "genomic_inverted_only"
    )

    s = d.loc[d.subset.eq(key)].copy()

    s["Study"] = s.study.map(STUDIES).fillna(s.study)

    s["Ranking method"] = s.method.map({
        **CHANNELS,
        "Prevalence": "Random reference",
        "Supervised": "Supervised ML",
    }).fillna(s.method)

    metric = st.radio(
        "Performance metric",
        ["AP", "ROC_AUC"],
        horizontal=True,
        key="benchmark_metric",
    )

    methods = [
        "ADPS",
        "Supervised",
        "case_expression_adps",
        "Prevalence",
        "EER",
        "dsRNAscan",
    ]

    selected = s.loc[s.method.isin(methods)]

    chart = selected.pivot(
        index="Study",
        columns="Ranking method",
        values=metric,
    )

    st.bar_chart(chart, stack=False)

    st.dataframe(
        s[[
            "Study", "Ranking method",
            "positive", "negative", "AP", "ROC_AUC"
        ]].round(3),
        hide_index=True,
        use_container_width=True,
    )

    st.caption(
        "The random AP reference is positive prevalence; "
        "chance ROC-AUC is 0.5."
    )

    st.warning(
        "Z22 is an expanded diagnostic benchmark: "
        "20 positives originate from five shared-TE-locus "
        "groups, and 17 of the 20 positives are "
        "genomically direct pairs. Its performance "
        "should not be interpreted as generalization "
        "to independent studies."
    )

    st.subheader("Additional benchmark analyses")

    with st.expander("Genomic-orientation counts"):
        st.dataframe(
            load_table(
                "01_benchmarks/ORIENTATION_COUNTS.csv"
            ),
            hide_index=True,
            use_container_width=True,
        )

    with st.expander("Paired statistical comparisons"):
        st.dataframe(
            load_table(
                "01_benchmarks/"
                "PAIRED_BOOTSTRAP_DIFFERENCES.csv"
            ),
            hide_index=True,
            use_container_width=True,
        )

    with st.expander("Search-radius sensitivity"):
        radius = load_table(
            "02_radius/TableS2_radius_metrics.csv"
        ).copy()

        radius["score_source"] = radius["score_source"].replace({
            "Prior Table S2 (rounded to 3 decimals)":
                "Radius-specific benchmark",
            "anchor_matched_user_verified_config":
                "Radius-specific benchmark",
        })

        radius = radius.rename(
            columns={"score_source": "Evaluation basis"}
        )

        st.dataframe(
            radius,
            hide_index=True,
            use_container_width=True,
        )

        st.caption(
            "All radii were evaluated without the "
            "dinucleotide-shuffle null model. "
            "FLAG-ZBP1, J2 and Z22 use radius-specific "
            "labeled cohorts. Table S2 reports their "
            "radius-specific counts and performance."
        )

    with st.expander("Cross-study transfer"):
        st.write(
            "Twenty directed transfers assess how "
            "models trained on one study perform "
            "when applied to another. Performance "
            "varies across assay contexts."
        )

        st.caption(
            "The Z22 GSE184962 transfer endpoint "
            "uses the broader case-enrichment cohort "
            "(28 positive and 38 negative), unlike "
            "the condition-specific principal benchmark."
        )

        st.dataframe(
            load_table(
                "03_transfer/"
                "TRANSFER_SEQUENCE_SCREEN_SUMMARY.csv"
            ).round(3),
            hide_index=True,
            use_container_width=True,
        )

        st.caption(
            "The sequence screen is a sensitivity "
            "analysis of fixed archived predictions, "
            "not a proof of absence of all homology."
        )

    with st.expander("Precision at top K"):
        st.dataframe(
            load_table(
                "01_benchmarks/TIE_AWARE_PRECISION_AT_K.csv"
            ),
            hide_index=True,
            use_container_width=True,
        )

    st.info(
        "ADPS is the label-independent ranking option. "
        "Supervised ranking requires assay-matched labels. "
        "A universal pretrained classifier is not "
        "recommended by the current transfer results."
    )

    st.markdown(
        f"[Download manuscript benchmark data]({DATA_DOI})"
    )


def documentation_page():
    st.header("Documentation")

    st.markdown(f"""
### Full dsRNASeeker workflow

The standalone pipeline accepts FASTQ or BAM input.

- [Installation and usage instructions]({GITHUB})
- [Archived software v1.0.4]({SOFTWARE_DOI})
- [Manuscript reproduction instructions]({GITHUB}/tree/main/reproduce_manuscript)
- [Benchmark and reproduction data]({DATA_DOI})

### Online analysis

1. Open the Run an analysis page.
2. Choose ADPS or assay-matched supervised ranking.
3. Try the bundled input example.
4. Upload a compatible candidate/evidence table.
5. Inspect and download the ranked results.

The bundled analysis examples are software test inputs.
The Candidate explorer instead displays actual archived
manuscript benchmark candidates.

ADPS is a prioritization score, not a calibrated probability
of physical dsRNA formation.
""")

    for filename in (
        "example_adps_input.tsv",
        "example_supervised_audit.tsv",
    ):
        st.download_button(
            f"Download {filename}",
            (EXAMPLES / filename).read_bytes(),
            file_name=filename,
            mime="text/tab-separated-values",
        )


def citation_page():
    st.header("Cite dsRNASeeker")

    st.subheader("Software")

    st.markdown(
        "dsRNASeeker v1.0.4 — Reproducibility "
        "and Cross-Study Benchmarking Update."
    )

    st.markdown(
        f"[DOI: 10.5281/zenodo.23407322]({SOFTWARE_DOI})"
    )

    st.subheader("Benchmark and reproduction data")

    st.write(
        "Levin, M. (2026). dsRNASeeker v1.0.4: "
        "Benchmarking and Manuscript Reproduction Data. "
        "Zenodo."
    )

    st.markdown(
        f"[DOI: 10.5281/zenodo.24296933]({DATA_DOI})"
    )

    st.markdown(f"[Source code and documentation]({GITHUB})")

    st.caption(
        "The November web-interface update is maintained "
        "in GitHub. The original software DOI identifies "
        "the frozen v1.0.4 software release."
    )


st.title("dsRNASeeker")

st.caption(
    "Double-stranded RNA candidate prioritization "
    "| Core software v1.0.4"
)

with st.sidebar:
    st.radio(
        "Navigation",
        PAGES,
        key="portal_page",
    )

    st.divider()

    st.markdown(f"[Source code]({GITHUB})")
    st.markdown(f"[Software archive]({SOFTWARE_DOI})")
    st.markdown(f"[Benchmark data]({DATA_DOI})")

    st.caption(
        "Computational evidence supports candidate "
        "prioritization, not definitive structural validation."
    )

routes = {
    "Home": home_page,
    "Candidate explorer": candidate_page,
    "Run an analysis": analysis_page,
    "Benchmarks": benchmark_page,
    "Documentation": documentation_page,
    "Citation": citation_page,
}

routes[st.session_state["portal_page"]]()
