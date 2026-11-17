#!/usr/bin/env python3
"""Plot completed FCCC results. No rescoring, re-fitting, or pipeline runs."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score

def main():
 p=argparse.ArgumentParser();p.add_argument('--support',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 s=a.support; m=s/'01_benchmarks'; r=s/'02_radius'
 studies=['FLAG_308489','J2_244103','Z22_184962']; titles=['FLAG-ZBP1 · 1 kb','J2 · 1 kb','Z22 · 5 kb (diagnostic)']; names=['FLAG 308489','J2 244103','Z22 184962']; cols=['#1f77b4','#ff7f0e','#238b69','#9467bd','#b33e38']; methods=['ADPS','Supervised','EER','dsRNAscan','case_expression_adps']; labels=['ADPS','Supervised ML','EER projected','dsRNAscan two-arm match','TE expression alone']; expected=[(75,52),(642,575),(20,38)]
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.titleweight':'bold','axes.labelsize':10,'legend.fontsize':8,'svg.fonttype':'none','pdf.fonttype':42})
 fig,ax=plt.subplots(3,2,figsize=(10.0,10.3),layout='constrained')
 for i,st in enumerate(studies):
  d=pd.read_csv(m/f'{st}_verified_labeled_scores.csv');y=d.label.astype(int);assert (int(y.sum()),len(y)-int(y.sum()))==expected[i]
  mm=pd.read_csv(m/'ALL_CHANNELS_AND_STRICT_SUBSET.csv').query('study==@st and subset=="all"').set_index('method')
  for method,label,col in zip(methods,labels,cols):
   ap=average_precision_score(y,d[method]);auc=roc_auc_score(y,d[method]);assert np.allclose([ap,auc],mm.loc[method,['AP','ROC_AUC']].astype(float),atol=1e-10,rtol=0)
   pr,re,_=precision_recall_curve(y,d[method]);fpr,tpr,_=roc_curve(y,d[method]);ax[i,0].plot(re,pr,color=col,lw=1.5,label=f'{label}: AP {ap:.3f}');ax[i,1].plot(fpr,tpr,color=col,lw=1.5,label=f'{label}: AUC {auc:.3f}')
  ax[i,0].axhline(y.mean(),color='.55',ls=':',lw=1.2,label=f'Prevalence: AP {y.mean():.3f}');ax[i,1].plot([0,1],[0,1],color='.55',ls=':',lw=1.2,label='Chance: AUC 0.500')
  for j,(xlab,ylab,kind) in enumerate([('Recall','Precision','Precision–recall'),('False-positive rate','True-positive rate','ROC')]):
   z=ax[i,j];z.set(title=f'{titles[i]} — {kind}',xlabel=xlab,ylabel=ylab,xlim=(0,1),ylim=(0,1.02),xticks=np.arange(0,1.01,.2),yticks=np.arange(0,1.01,.2));z.grid(alpha=.15);z.legend(loc='lower right',frameon=True,framealpha=.93,edgecolor='none',fontsize=7.4)
 for e in ['png','pdf','svg']: fig.savefig(a.out/f'Figure2_FINAL.{e}',dpi=350)
 plt.close(fig)
 print('PASS: Figure 2: six lower-right legends; candidate scores verified.')
if __name__=='__main__':main()
