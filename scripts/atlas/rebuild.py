#!/usr/bin/env python3
"""Rebuild atlas quantities from all 641 Salmon archives using final frozen curation.

Raw TPM is summed across transcripts per gene; BioSample medians precede study
medians and the median of study medians. No cross-study inferential test is made.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parents[2];B=P/'source_addenda/B_atlas_salmon';O=P/'results/atlas_rebuilt';O.mkdir(parents=True,exist_ok=True)
a=argparse.ArgumentParser();a.add_argument('--matrix',action='store_true',help='Also write the complete gene TPM matrix (~45 MB compressed)');args=a.parse_args()
FOCAL='Vitvi05_01chr01g05750'
REFS={'UBQ2':'Vitvi05_01chr19g07510','UBQ':'Vitvi05_01chr19g07460','PP2A-A':'Vitvi05_01chr17g05170','GAPDH':'Vitvi05_01chr17g15130','EF1a':'Vitvi05_01chr06g04820','SAND':'Vitvi05_01chr06g04330','ACT1':'Vitvi05_01chr04g23050'}
per=pd.read_csv(P/'inputs/atlas/focal_atlas_per_run.tsv',sep='\t');per=per.set_index('run_accession');old=per.copy();mat={}
for i,f in enumerate(sorted((B/'quant').glob('*/quant.sf.gz'))):
 run=f.parent.name;q=pd.read_csv(f,sep='\t',usecols=['Name','TPM','NumReads']);q['gene']=q.Name.str.replace(r'_t\d+$','',regex=True)
 t=q.groupby('gene',sort=True).TPM.sum();m=json.loads((f.parent/'meta_info.json').read_text())
 per.loc[run,'focal_tpm']=t.loc[FOCAL];per.loc[run,'focal_reads']=q.loc[q.gene==FOCAL,'NumReads'].sum()
 per.loc[run,'n_processed']=m['num_processed'];per.loc[run,'pct_mapped']=m['percent_mapped']
 for k,g in REFS.items():per.loc[run,k+'_tpm']=t.loc[g]
 if args.matrix:mat[run]=t
 if (i+1)%100==0:print('Quantifications read:',i+1,flush=True)
assert len(per)==641
maxerr=float((per.focal_tpm-old.focal_tpm).abs().max());assert maxerr<1e-6,maxerr
assert float((per.focal_reads-old.focal_reads).abs().max())<1e-6
per['focal_cpm']=per.focal_reads/per.n_processed*1e6
per['log2_focal_over_UBQ2']=np.log2((per.focal_tpm+.5)/(per.UBQ2_tpm+.5))
per['qc_ok']=(per.n_processed>=5e6)&(per.pct_mapped>=40)
per.reset_index().to_csv(O/'focal_atlas_per_run.tsv',sep='\t',index=False)
selected=per[per.qc_ok & per.in_selection & per.tissue_class.notna()].reset_index()
vals=['focal_tpm','focal_cpm','log2_focal_over_UBQ2']+[k+'_tpm' for k in REFS]
def summarize(rs,prefix):
 s=rs.groupby(['tissue_class','study_accession','sample_accession']).agg(n_runs=('run_accession','size'),cultivar=('cultivar','first'),class_basis=('class_basis','first'),submitter_class=('submitter_class','first'),**{v:(v,'median') for v in vals}).reset_index()
 agg={'n_samples':('sample_accession','size'),'n_runs':('n_runs','sum'),'median_focal_tpm':('focal_tpm','median'),'median_focal_cpm':('focal_cpm','median'),'class_basis':('class_basis','first'),'submitter_class':('submitter_class','first'),**{f'median_{v}':(v,'median') for v in vals if v not in ['focal_tpm','focal_cpm']}}
 b=s.groupby(['tissue_class','study_accession']).agg(**agg).reset_index()
 t=b.groupby('tissue_class').agg(studies=('study_accession','nunique'),samples=('n_samples','sum'),runs=('n_runs','sum'),median_study_focal_tpm=('median_focal_tpm','median'),median_study_focal_cpm=('median_focal_cpm','median')).reset_index()
 for name,d in [('per_sample',s),('by_study',b),('by_tissue',t)]:d.to_csv(O/f'{prefix}focal_atlas_{name}.tsv',sep='\t',index=False)
 return s,b,t
s,_,_=summarize(selected,'');primary=selected[selected.tissue_class!='unlabelled calli (inferred embryogenic)'];sp,_,_=summarize(primary,'primary_')
assert (len(primary),len(sp),primary.study_accession.nunique())==(624,591,48)
assert (len(selected),len(s),selected.study_accession.nunique())==(630,597,49)
saved=pd.read_csv(P/'inputs/atlas/focal_atlas_per_sample.tsv',sep='\t')
cmp=s.merge(saved,on=['tissue_class','study_accession','sample_accession'],suffixes=('_new','_saved'),validate='one_to_one');assert len(cmp)==597
assert np.allclose(cmp.focal_tpm_new,cmp.focal_tpm_saved,atol=1e-6,rtol=1e-12)
if args.matrix:pd.DataFrame(mat).to_csv(O/'gene_tpm.tsv.gz',sep='\t',float_format='%.3f')
(O/'rebuild_checks.json').write_text(json.dumps({'quantifications':641,'primary_runs':624,'primary_BioSamples':591,'primary_studies':48,'full_runs':630,'full_BioSamples':597,'full_studies':49,'max_focal_TPM_discrepancy':maxerr,'BioSample_medians_match':True,'matrix_written':args.matrix},indent=2)+'\n')
print('PASS: all 641 focal quantities and 597 final BioSample medians match frozen atlas inputs.')
