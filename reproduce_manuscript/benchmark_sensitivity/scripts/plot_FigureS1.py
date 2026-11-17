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
 p=argparse.ArgumentParser()
 p.add_argument('--support',type=Path,required=True)
 p.add_argument('--out',type=Path,required=True)
 a=p.parse_args()
 a.out.mkdir(parents=True,exist_ok=True)

 r=a.support/'02_radius'
 d=pd.read_csv(r/'TableS2_radius_metrics.csv')
 rho=pd.read_csv(r/'TableS2_shared_rank_stability.csv')

 names=['FLAG 308489','J2 244103','Z22 184962']

 plt.rcParams.update({
 'font.family':'DejaVu Sans',
 'font.size':10,
 'axes.titlesize':11,
 'axes.titleweight':'bold',
 'axes.labelsize':10,
 'legend.fontsize':8,
 'svg.fonttype':'none',
 'pdf.fonttype':42})
 fig,ax=plt.subplots(2,2,figsize=(10,8.5));fig.subplots_adjust(left=.09,right=.98,bottom=.08,top=.90,hspace=.30,wspace=.19)
 for st,n,col,mark in zip(['FLAG-ZBP1','J2','Z22'],names,['#1f77b4','#ff7f0e','#2ca02c'],['o','s','^']):
  q=d[d.study.eq(st)];x=q.radius_bp/1000;ax[0,0].plot(x,q.candidate_n,marker=mark,color=col,lw=1.8,ms=7,label=n)
  for z,field in [(ax[0,1],'AP'),(ax[1,0],'ROC_AUC')]:
   if st=='Z22':
    z.scatter(x.iloc[:2],q[field].iloc[:2],marker=mark,s=100,facecolors='none',edgecolors=col,linewidths=1.8);z.scatter(x.iloc[2:],q[field].iloc[2:],marker=mark,s=100,color=col)
   else:z.scatter(x,q[field],marker=mark,s=100,color=col)
 ax[0,0].set(yscale='log',ylabel='Candidate pairs (log scale)',title='A  Candidate-pair enumeration');ax[0,1].set(ylabel='AP',title='B  Average precision',ylim=(-.03,1.02));ax[1,0].set(ylabel='ROC-AUC',title='C  ROC-AUC',ylim=(-.03,1.02));ax[1,0].axhline(.5,color='.55',ls='--',lw=1)
 for z in [ax[0,0],ax[0,1],ax[1,0]]:z.set(xlabel='Maximum TE-pair distance (kb)',xticks=[1,2,5]);z.title.set_ha('left');z.title.set_position((0,1.0))
 q=rho[(rho.radius_1==1000)&(rho.radius_2==5000)].set_index('study').loc[['FLAG-ZBP1','J2','Z22']];bars=ax[1,1].bar(names,q.spearman,color='#1f77b4',width=.75);ax[1,1].set(ylim=(0,1.1),yticks=np.arange(0,1.01,.2),ylabel='Spearman ρ',title='D  Shared-pair rank stability: 1 vs 5 kb');ax[1,1].bar_label(bars,fmt='%.4f',padding=3);ax[1,1].title.set_ha('left');ax[1,1].title.set_position((0,1.0));ax[1,1].tick_params(axis='x',labelsize=9)
 fig.legend(*ax[0,0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.5,.98),fontsize=11)
 for e in ['png','pdf','svg']:fig.savefig(a.out/f'FigureS1_FINAL.{e}',dpi=350)
 plt.close(fig)

 print('PASS: Figure S1 recreated from archived radius results.')

if __name__=='__main__':
 main()
