from __future__ import annotations
from io import BytesIO
import pandas as pd
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.responses import JSONResponse
from .scoring import calculate_adps
from .supervised import same_study_rerank
from .registry import load_json

app=FastAPI(title="dsRNASeeker Web API",version="0.1.0")

def read_table(raw:bytes,name:str)->pd.DataFrame:
    sep="\t" if name.lower().endswith((".tsv",".txt")) else ","
    return pd.read_csv(BytesIO(raw),sep=sep,low_memory=False)

@app.get("/health")
def health(): return {"status":"ok"}
@app.get("/benchmarks")
def benchmarks(): return load_json("studies.json")
@app.get("/models")
def models(): return load_json("models.json")

@app.post("/score/adps")
async def score_adps(file:UploadFile=File(...),case:str=Form(...),control:str=Form(...),annotation_policy:str=Form("conservative")):
    try:
        df=read_table(await file.read(),file.filename or "input.tsv")
        out=calculate_adps(df,case,control,annotation_policy)
        return JSONResponse({"n":len(out),"columns":list(out.columns),"rows":out.head(500).where(pd.notna(out),None).to_dict("records")})
    except Exception as e: raise HTTPException(status_code=400,detail=str(e))

@app.post("/rerank/same-study")
async def rerank(file:UploadFile=File(...),outer_folds:int=Form(5),seed:int=Form(20260925)):
    try:
        df=read_table(await file.read(),file.filename or "audit.tsv")
        ranked,folds=same_study_rerank(df,outer_folds,seed)
        return JSONResponse({"n":len(ranked),"ranking":ranked.head(500).where(pd.notna(ranked),None).to_dict("records"),"folds":folds.where(pd.notna(folds),None).to_dict("records")})
    except Exception as e: raise HTTPException(status_code=400,detail=str(e))
