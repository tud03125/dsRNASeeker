#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from sklearn.metrics import average_precision_score, roc_auc_score

SPECS = [
    ("GSE308489__w1000_standardized", "GSE308488/9\nFLAG-ZBP1 · 1 kb"),
    ("GSE184962_4__w5000_expanded", "GSE184962/4\nZ22 · 5 kb"),
    ("GSE244103__w1000_standardized", "GSE244103\nJ2 · 1 kb"),
]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--outdir',required=True)
    a=ap.parse_args(); root=Path(a.root); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    metrics=[]; retention=[]; audits=[]
    for ds,label in SPECS:
        audit_p=root/'supervised_input_expanded'/f'{ds}.audit.tsv'
        map_p=root/'mapping_confidence'/ds/'candidate_mapping_confidence.tsv'
        audit=pd.read_csv(audit_p,sep='\t',low_memory=False)
        mc=pd.read_csv(map_p,sep='\t',low_memory=False)
        q=audit[['pair_id','label','adaptive_priority_score']].merge(
            mc[['pair_id','label','high_mapping_confidence','ADPS']], on='pair_id', how='inner', suffixes=('_audit','_map'), validate='one_to_one')
        if len(q)!=len(audit) or len(q)!=len(mc):
            raise RuntimeError(f'{ds}: pair-id join not one-to-one complete: audit={len(audit)} map={len(mc)} join={len(q)}')
        la=pd.to_numeric(q.label_audit,errors='coerce'); lm=pd.to_numeric(q.label_map,errors='coerce')
        labok=la.notna() & lm.notna()
        if not np.array_equal(la[labok].to_numpy(),lm[labok].to_numpy()):
            raise RuntimeError(f'{ds}: label mismatch between audit and mapping tables')
        score=pd.to_numeric(q.adaptive_priority_score,errors='coerce')
        old=pd.to_numeric(q.ADPS,errors='coerce')
        both=score.notna() & old.notna()
        maxdiff=float((score[both]-old[both]).abs().max()) if both.any() else np.nan
        audits.append({'dataset':ds,'n_pairs':len(q),'max_abs_difference_canonical_vs_mapping_ADPS':maxdiff,'n_difference_gt_1e-12':int(((score[both]-old[both]).abs()>1e-12).sum())})
        q['label']=la
        q['canonical_adps']=score
        q['high_mapping_confidence']=q.high_mapping_confidence.astype(str).str.lower().isin(['1','true','t','yes','y'])
        for class_name,mask in [
            ('positive',q.label.eq(1)),('confident_negative',q.label.eq(0)),('unlabeled',q.label.isna()),('all_candidates',pd.Series(True,index=q.index))]:
            n0=int(mask.sum()); nh=int((mask & q.high_mapping_confidence).sum())
            retention.append({'dataset':ds,'display_dataset':label.replace('\n',' '),'label_class':class_name,'n_before':n0,'n_high_mapping_confidence':nh,'retention_fraction':nh/n0 if n0 else np.nan,'removed':n0-nh})
        labeled=q[q.label.isin([0,1]) & q.canonical_adps.notna()].copy()
        for subset,keep in [('all_labeled',pd.Series(True,index=labeled.index)),('high_mapping_confidence',labeled.high_mapping_confidence)]:
            d=labeled.loc[keep]
            y=d.label.astype(int); s=d.canonical_adps.astype(float)
            metrics.append({'dataset':ds,'display_dataset':label.replace('\n',' '),'subset':subset,'n':len(d),'n_positive':int(y.sum()),'n_negative':int((1-y).sum()),'average_precision':average_precision_score(y,s),'roc_auc':roc_auc_score(y,s)})
    pd.DataFrame(audits).to_csv(out/'mapping_confidence_canonical_score_audit.tsv',sep='\t',index=False)
    pd.DataFrame(retention).to_csv(out/'mapping_confidence_retention_by_label.tsv',sep='\t',index=False)
    met=pd.DataFrame(metrics); met.to_csv(out/'mapping_confidence_canonical_ADPS_metrics.csv',index=False)

    fig,axes=plt.subplots(1,2,figsize=(12,5.6))
    x=np.arange(len(SPECS)); width=0.34
    for ax,metric,ylabel in [(axes[0],'average_precision','Average precision'),(axes[1],'roc_auc','ROC-AUC')]:
        for j,(subset,legend) in enumerate([('all_labeled','All labeled'),('high_mapping_confidence','High mapping confidence')]):
            s=met[met.subset.eq(subset)].set_index('dataset').loc[[z[0] for z in SPECS]]
            ax.bar(x+(j-0.5)*width,s[metric].to_numpy(),width,label=legend)
        ax.set_xticks(x,[z[1] for z in SPECS],rotation=0)
        ax.set_ylim(0,1.02); ax.set_ylabel(ylabel); ax.grid(axis='y',alpha=.2)
        if metric=='roc_auc': ax.axhline(.5,linestyle='--',linewidth=1,alpha=.65)
    axes[0].set_title('A  Ranking performance',loc='left')
    axes[1].set_title('B  Discrimination performance',loc='left')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.suptitle('ADPS benchmark performance after stringent mapping-confidence filtering',y=.98,fontsize=15)
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.925),ncol=2,frameon=False)
    fig.subplots_adjust(left=.08,right=.985,top=.82,bottom=.17,wspace=.18)
    for ext in ['png','pdf','svg']:
        fig.savefig(out/f'FigureS_mapping_confidence_canonical.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
    print(pd.DataFrame(audits).to_string(index=False))
    print(met.to_string(index=False))

if __name__=='__main__': main()
