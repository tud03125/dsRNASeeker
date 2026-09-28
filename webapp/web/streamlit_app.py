from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import streamlit as st
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from app.scoring import calculate_adps
from app.supervised import same_study_rerank
from app.registry import load_json

st.set_page_config(page_title="dsRNASeeker Web Explorer",layout="wide")
st.title("dsRNASeeker Web Explorer")
st.caption("Prototype: candidate-table ADPS scoring, assay-matched supervised re-ranking, and frozen benchmark context. No universal pretrained classifier is applied by default.")

t1,t2,t3=st.tabs(["ADPS scoring","Same-study supervised","Benchmark registry"])
with t1:
    f=st.file_uploader("Upload dsRNASeeker candidate/evidence CSV or TSV",key="adps")
    case=st.text_input("Case label",value="CASE"); control=st.text_input("Control label",value="CONTROL")
    policy=st.selectbox("Annotation policy",["conservative","zrna_permissive"])
    if f and st.button("Calculate ADPS"):
        df=pd.read_csv(f,sep="\t" if f.name.endswith((".tsv",".txt")) else ",",low_memory=False)
        try:
            out=calculate_adps(df,case,control,policy); st.dataframe(out.head(200),use_container_width=True)
            st.download_button("Download ranked TSV",out.to_csv(sep="\t",index=False).encode(),"dsRNASeeker_ranked.tsv")
        except Exception as e: st.error(str(e))
with t2:
    st.info("Use assay-matched labels. Current cross-study benchmarks do not support a universal pretrained classifier.")
    f=st.file_uploader("Upload labeled audit TSV/CSV",key="sup")
    if f and st.button("Run grouped same-study re-ranking"):
        df=pd.read_csv(f,sep="\t" if f.name.endswith((".tsv",".txt")) else ",",low_memory=False)
        try:
            ranked,folds=same_study_rerank(df); st.dataframe(folds,use_container_width=True); st.dataframe(ranked.head(200),use_container_width=True)
            st.download_button("Download supervised ranking TSV",ranked.to_csv(sep="\t",index=False).encode(),"dsRNASeeker_supervised_ranked.tsv")
        except Exception as e: st.error(str(e))
with t3:
    st.json(load_json("studies.json")); st.json(load_json("models.json"))
