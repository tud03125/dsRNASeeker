#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

SPECS=[
('GSE308489__w1000_standardized','w1000_standardized','GSE308488/9\nFLAG-ZBP1 · 1 kb'),
('GSE184962_4__w5000_expanded','w5000_expanded','GSE184962/4\nZ22 · 5 kb'),
('GSE244103__w1000_standardized','w1000_standardized','GSE244103\nJ2 · 1 kb'),
]
SERIES=[
('same_study_fixed_baseline','integrated_adps','ADPS','o'),
('same_study_nested_grouped_cv','fixed_l2_raw7','Same-study supervised','s'),
('leave_one_study_family_out','loso_regularized_prespecified','Held-out-family supervised','^'),
]

def one(df,ds,var,ev,model):
 q=df[(df.dataset==ds)&(df.analysis_variant==var)&(df.evaluation==ev)&(df.model==model)]
 if len(q)!=1: raise RuntimeError(f'expected one row for {ds} {var} {ev} {model}, got {len(q)}')
 return q.iloc[0]

def main():
 p=argparse.ArgumentParser(); p.add_argument('--table',required=True); p.add_argument('--outdir',required=True)
 a=p.parse_args(); d=pd.read_csv(a.table); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
 fig,axes=plt.subplots(1,2,figsize=(13.2,6.6)); x=np.arange(len(SPECS)); offsets=[-.17,0,.17]
 for ax,metric,lo,hi,ylabel,panel in [
  (axes[0],'average_precision','ap_ci_low','ap_ci_high','Average precision','A  Ranking performance'),
  (axes[1],'roc_auc','auc_ci_low','auc_ci_high','ROC-AUC','B  Discrimination performance')]:
  for j,(ev,model,name,marker) in enumerate(SERIES):
   vals=[]; lows=[]; highs=[]
   for ds,var,label in SPECS:
    r=one(d,ds,var,ev,model); v=float(r[metric]); l=float(r[lo]); h=float(r[hi])
    vals.append(v); lows.append(max(0,v-l)); highs.append(max(0,h-v))
   xx=x+offsets[j]
   ax.errorbar(xx,vals,yerr=[lows,highs],fmt=marker,capsize=4,markersize=7,linewidth=1.5,label=name)
  ax.set_ylim(0,1.02); ax.set_xticks(x,[z[2] for z in SPECS],fontsize=10); ax.set_ylabel(ylabel,fontsize=12); ax.grid(axis='y',alpha=.2); ax.set_title(panel,loc='left',fontsize=14)
  if metric=='roc_auc': ax.axhline(.5,linestyle='--',linewidth=1.1,alpha=.65)
 handles,labels=axes[0].get_legend_handles_labels()
 fig.suptitle('Same-study adaptation does not reliably transfer across assay families',fontsize=16,y=.975)
 fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.915),ncol=3,frameon=False,fontsize=10.5)
 fig.text(.5,.025,'Z22 uses the 5-kb labeled representative (20 positive / 38 confident-negative); the 1-kb subset (5 / 3) does not meet the ≥10 / ≥10 same-study eligibility rule.',ha='center',fontsize=9)
 fig.subplots_adjust(left=.08,right=.985,top=.80,bottom=.17,wspace=.18)
 for ext in ['png','pdf','svg']:
  fig.savefig(out/f'Figure4_supervised_generalization_fixed.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
 print('[OK]',out/'Figure4_supervised_generalization_fixed.pdf')
if __name__=='__main__': main()
