#!/usr/bin/env python3
"""Paired group bootstrap of saved scores; full and genomic-inverted-only subsets.
No feature recomputation or supervised re-training; supervised scores are pooled OOF.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.metrics import average_precision_score,roc_auc_score,precision_recall_curve,roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
F=['orientation_adps','annotation_adps','case_expression_adps','energy_adps','interface_adps','case_editing_adps','RI_adps']
NAMES=dict(zip(F,['Orientation','Annotation','Case TE expression','Duplex energy','Pairing interface','A-to-I editing','Intron retention']))
METHODS=['ADPS','Supervised','EER','dsRNAscan','dsRNAscan_stability','dsRNAscan_probing']+F+['Prevalence']
PRIMARY=['EER','dsRNAscan']; SECONDARY=['case_expression_adps']
STUDIES=['FLAG_308489','J2_244103','Z22_184962']
class UF:
    def __init__(self):self.p={}
    def find(self,x):
        self.p.setdefault(x,x);r=x
        while self.p[r]!=r:r=self.p[r]
        while self.p[x]!=x:y=self.p[x];self.p[x]=r;x=y
        return r
    def union(self,a,b):self.p[self.find(a)]=self.find(b)
def groups(d):
    u=UF()
    for r in d.itertuples(index=False):
        a,b='arm:'+r.A_TE_id,'arm:'+r.B_TE_id;g='group:'+str(r.group_id)
        u.union(a,b);u.union(a,g)
    return d.A_TE_id.map(lambda a:u.find('arm:'+a))
def metrics(y,s):return np.array([average_precision_score(y,s),roc_auc_score(y,s)])
def wmetrics(y,s,w):
    # Vectorized weighted AP and AUC with exact tied-score handling.
    order=np.argsort(-s,kind='stable'); y=y[order];s=s[order];w=w[:,order]
    starts=np.r_[0,np.flatnonzero(s[1:]!=s[:-1])+1]
    p=np.add.reduceat(w*y,starts,axis=1);n=np.add.reduceat(w*(1-y),starts,axis=1)
    cp=p.cumsum(axis=1);cn=n.cumsum(axis=1);P=cp[:,-1];N=cn[:,-1]
    pr=np.divide(cp,cp+cn,out=np.zeros_like(cp,dtype=float),where=(cp+cn)>0)
    ap=(p*pr).sum(axis=1)/P
    auc=(p*(N[:,None]-cn+.5*n)).sum(axis=1)/(P*N)
    return np.column_stack([ap,auc])
def resample(d,cols,gcol,B,seed):
    y=d.label.to_numpy(dtype=int); g,ix=np.unique(d[gcol].astype(str),return_inverse=True);rng=np.random.default_rng(seed)
    W=rng.multinomial(len(g),np.full(len(g),1/len(g)),size=B)[:,ix]; valid=(W@y>0)&(W@(1-y)>0);W=W[valid];assert len(W)>0
    out={}
    for c in cols:
        s=d[c].to_numpy(float);out[c]=np.concatenate([wmetrics(y,s,W[i:i+256]) for i in range(0,len(W),256)])
        # Check first resample against independent sklearn implementation.
        assert np.allclose(out[c][0],[average_precision_score(y,s,sample_weight=W[0]),roc_auc_score(y,s,sample_weight=W[0])],atol=1e-12)
    return len(g),len(W),out

def holm(p):
    p=np.asarray(p,float);o=np.argsort(p);v=np.maximum.accumulate((len(p)-np.arange(len(p)))*p[o]);z=np.empty_like(p);z[o]=np.minimum(v,1);return z

def main():
    a=argparse.ArgumentParser();a.add_argument('--inputs',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--bootstrap',type=int,default=5000);a.add_argument('--seed',type=int,default=20261114);q=a.parse_args();q.out.mkdir(parents=True,exist_ok=True);assert q.bootstrap>=100
    rows=[];diffs=[];counts=[];groupsrows=[];data={};pts=[]; B=q.bootstrap
    for i,study in enumerate(STUDIES):
        full=pd.read_csv(q.inputs/f'{study}.csv.gz');assert full.pair_id.is_unique;full['connected_group']=groups(full);full['Prevalence']=0.0
        for orient,g in full.groupby('genomic_orientation'):
            counts.append(dict(study=study,orientation=orient,candidates=len(g),positives=int(g.label.eq(1).sum()),negatives=int(g.label.eq(0).sum()),unlabeled=int(g.label.isna().sum())))
        all_lab=full[full.label.isin([0,1])].copy();data[study]=all_lab
        for subset,d in [('all',all_lab),('genomic_inverted_only',all_lab[all_lab.genomic_orientation.eq('inverted')])]:
            d=d.copy();y=d.label.to_numpy(int);assert len(set(y))==2;assert np.isfinite(d[METHODS].to_numpy(float)).all()
            for c in METHODS:
                ap,auc=metrics(y,d[c]); rows.append(dict(study=study,subset=subset,method=c,n=len(d),positive=int(y.sum()),negative=int(len(y)-y.sum()),AP=ap,ROC_AUC=auc,AP_CI_low=np.nan,AP_CI_high=np.nan,AUC_CI_low=np.nan,AUC_CI_high=np.nan))
            # Existing scores only; strict subset does not refit weights or models.
            # Tie-aware precision at k on labeled cohorts, not all discovered candidates.
            for c in ['ADPS','Supervised','EER','dsRNAscan','case_expression_adps']:
                s=d[c].to_numpy(float)
                for k in [10,20]:
                    if k>len(d):continue
                    cut=np.sort(s)[-k];hi=s>cut;tie=s==cut;slots=k-int(hi.sum());p0=int(y[hi].sum());p1=int(y[tie].sum());t=int(tie.sum())
                    pts.append(dict(study=study,subset=subset,method=c,k=k,expected=(p0+slots*p1/t)/k,lower=(p0+max(0,slots-(t-p1)))/k,upper=(p0+min(slots,p1))/k))
            if subset!='all':continue
            for j,gcol in enumerate(['connected_group','group_id']):
                G,V,boot=resample(d,METHODS,gcol,B,q.seed+i*10+j);gp=d.loc[d.label.eq(1),gcol].nunique();gn=d.loc[d.label.eq(0),gcol].nunique()
                groupsrows.append(dict(study=study,scheme=gcol,groups=G,positive_groups=gp,negative_groups=gn,requested=B,valid=V,largest_group_n=int(d[gcol].value_counts().max())))
                if gcol=='connected_group':
                    for row in rows:
                        if row['study']==study and row['subset']=='all':
                            ci=np.quantile(boot[row['method']],[.025,.975],axis=0);row.update(AP_CI_low=ci[0,0],AP_CI_high=ci[1,0],AUC_CI_low=ci[0,1],AUC_CI_high=ci[1,1])
                for c in PRIMARY+SECONDARY:
                    delta=metrics(y,d.ADPS)-metrics(y,d[c]);b=boot['ADPS']-boot[c]
                    for ix,m in enumerate(['AP','ROC_AUC']):
                        lo,hi=np.quantile(b[:,ix],[.025,.975]);pv=(1+np.sum(np.abs(b[:,ix]-delta[ix])>=abs(delta[ix])))/(V+1)
                        diffs.append(dict(study=study,group_scheme=gcol,comparator=c,metric=m,delta=delta[ix],ci_low=lo,ci_high=hi,p_centered_bootstrap_approx=pv,valid_resamples=V,positive_groups=gp))
        all_lab.to_csv(q.out/f'{study}_verified_labeled_scores.csv',index=False)
    df=pd.DataFrame(rows);df.to_csv(q.out/'ALL_CHANNELS_AND_STRICT_SUBSET.csv',index=False);pd.DataFrame(counts).to_csv(q.out/'ORIENTATION_COUNTS.csv',index=False);pd.DataFrame(groupsrows).to_csv(q.out/'BOOTSTRAP_GROUPS.csv',index=False);d=pd.DataFrame(diffs);d['holm_p_12_primary']=np.nan
    for scheme in d.group_scheme.unique():
        mask=d.group_scheme.eq(scheme)&d.comparator.isin(PRIMARY);d.loc[mask,'holm_p_12_primary']=holm(d.loc[mask,'p_centered_bootstrap_approx'])
    d.to_csv(q.out/'PAIRED_BOOTSTRAP_DIFFERENCES.csv',index=False);pd.DataFrame(pts).to_csv(q.out/'TIE_AWARE_PRECISION_AT_K.csv',index=False)
    # Main figure: retain four published curves, add expression baseline and chance lines.
    # Other six channels in separate figure/table to avoid unreadable main plots.
    colors={'ADPS':'#2468A0','Supervised':'#E69F00','EER':'#157A60','dsRNAscan':'#8055A0','case_expression_adps':'#B33E38'}
    fig,ax=plt.subplots(3,2,figsize=(9,11),layout='constrained')
    for i,study in enumerate(STUDIES):
        d=data[study];y=d.label.to_numpy(int)
        for method in colors:
            s=d[method];ap,auc=metrics(y,s);p,r,_=precision_recall_curve(y,s);f,t,_=roc_curve(y,s);name=NAMES.get(method,method)
            ax[i,0].step(r,p,where='post',label=f'{name}: {ap:.3f}',color=colors[method],lw=1.5);ax[i,1].plot(f,t,label=f'{name}: {auc:.3f}',color=colors[method],lw=1.5)
        ax[i,0].axhline(y.mean(),color='0.45',ls=':',label=f'Prevalence: {y.mean():.3f}');ax[i,1].plot([0,1],[0,1],color='0.45',ls=':',label='Chance: 0.500')
        for j in [0,1]:
            ax[i,j].set(xlim=(0,1),ylim=(0,1.02),xlabel='Recall' if j==0 else 'False-positive rate',ylabel='Precision' if j==0 else 'True-positive rate',title=f'{study}'+(' (diagnostic)' if i==2 else '')+(' · PR' if j==0 else ' · ROC'));ax[i,j].legend(fontsize=7,loc='lower left' if j==0 else 'lower right')
    for ext in ['pdf','svg','png']:fig.savefig(q.out/f'Figure2_revised.{ext}',dpi=250)
    plt.close(fig)
    # All single channels together (supplement): same metric definitions, no fitted direction.
    fig,ax=plt.subplots(3,2,figsize=(9,9),layout='constrained');labels=[NAMES.get(c,c) for c in ['ADPS']+F]
    for i,study in enumerate(STUDIES):
        t=df.query('study==@study and subset=="all"').set_index('method').loc[['ADPS']+F]
        for j,m in enumerate(['AP','ROC_AUC']):
            ax[i,j].barh(labels,t[m],color=['#2468A0']+['#8CA5B5']*7);ax[i,j].invert_yaxis();ax[i,j].axvline(data[study].label.mean() if j==0 else .5,ls=':',color='0.35');ax[i,j].set(xlim=(0,1),title=f'{study} · {m}');ax[i,j].tick_params(labelsize=8)
    fig.savefig(q.out/'All_single_channels.pdf');fig.savefig(q.out/'All_single_channels.svg');plt.close(fig)
    meta=dict(seed=q.seed,requested_bootstrap=B,groups_primary='Union of benchmark-assigned groups and shared TE arms across complete candidate set, before restricting to labeled pairs',group_sensitivity='benchmark group_id alone',inferential_scope='paired fixed-score cluster bootstrap, pointwise percentile 95% CI; approximate centered-bootstrap two-sided p, Holm over 12 ADPS vs EER/two-arm dsRNAscan endpoints within each group scheme; expression comparison exploratory, no Holm value',strict_subset='restricted evaluation of existing scores and OOF predictions; no refitting or recalculation of ADPS weights',baseline='constant score; AP positive prevalence, AUC 0.5; not simulated random permutation',input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in q.inputs.glob('*.csv.gz')})
    (q.out/'ANALYSIS_METADATA.json').write_text(json.dumps(meta,indent=2)+'\n');print(pd.DataFrame(groupsrows).to_string(index=False));print('Saved:',q.out)
if __name__=='__main__':main()
