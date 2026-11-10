import csv,gzip,json
from pathlib import Path
import numpy as np
import os
root=Path('.');out=root/'results/paralog';out.mkdir(parents=True,exist_ok=True)
co=root/os.environ.get('VVECPG1_RESULTS','results')/'coexpression'
focal='Vitvi05_01chr01g05750';partner='Vitvi05_01chr01g05760';genes={focal,partner}
def rows(path):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt') as f:
  r=csv.DictReader(f,delimiter='\t');key=r.fieldnames[0];return {d[key]:d for d in r if d[key] in genes}
raw=rows('inputs/counts/pn_v51.primary.counts.tsv');norm=rows('inputs/counts/pn_v51.normalized_counts_deseq2.tsv');atlas=rows('inputs/atlas/gene_tpm.tsv.gz')
with open('inputs/metadata/samples.tsv') as f:meta={d['sample_uid']:d for d in csv.DictReader(f,delimiter='\t')}
with open(co/'discovery_samples.tsv') as f:dis=[d['sample_uid'] for d in csv.DictReader(f,delimiter='\t')]
with open(co/'validation_samples.tsv') as f:val=[d['run_accession'] for d in csv.DictReader(f,delimiter='\t')]
sets={'all_126_libraries':list(meta),'early_90_treated':dis,'24h_30_treated':[s for s,d in meta.items() if d['time_hours']=='24' and d['is_control']=='False'],'7d_30_treated':[s for s,d in meta.items() if d['time_hours']=='168' and d['is_control']=='False'],'controls_6_preparations':[s for s,d in meta.items() if d['is_control']=='True']}
summary=[]
def summarize(context,unit,data,ss):
 x=np.array([float(data[focal][s]) for s in ss]);y=np.array([float(data[partner][s]) for s in ss]);r=x/y
 d=dict(context=context,unit=unit,n=len(ss),focal_median=float(np.median(x)),partner_median=float(np.median(y)),focal_min=float(x.min()),focal_max=float(x.max()),partner_min=float(y.min()),partner_max=float(y.max()),paired_ratio_median=float(np.median(r)),paired_ratio_min=float(r.min()),paired_ratio_max=float(r.max()),paired_ratio_q25=float(np.quantile(r,.25)),paired_ratio_q75=float(np.quantile(r,.75)),all_focal_greater=bool(np.all(x>y)))
 summary.append(d)
 for s,a,b in zip(ss,x,y):
  pass
for context,ss in sets.items():
 if not ss:continue
 for unit,data in [('raw_UMI_deduplicated_gene_counts',raw),('DESeq2_size_factor_normalized_counts',norm)]:summarize(context,unit,data,ss)
summarize('independent_27_embryogenic_samples','Salmon_gene_TPM',atlas,val)
with open(out/'expression_summary.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
with open(out/'hormone_pair_counts.tsv','w') as f:
 fields=['sample_uid','group','time_hours','batch','is_control','focal_raw','partner_raw','focal_normalized','partner_normalized','paired_count_ratio'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for s,d in meta.items():
  w.writerow(dict(sample_uid=s,group=d['group'],time_hours=d['time_hours'],batch=d['batch'],is_control=d['is_control'],focal_raw=raw[focal][s],partner_raw=raw[partner][s],focal_normalized=norm[focal][s],partner_normalized=norm[partner][s],paired_count_ratio=float(raw[focal][s])/float(raw[partner][s])))
for d in summary:
 if d['unit']!='raw_UMI_deduplicated_gene_counts':print(json.dumps(d))
# Explicit current candidate result, without approximate matching by other gene suffixes.
with open(co/'fixed_candidates_validation.tsv') as f:
 d=next(d for d in csv.DictReader(f,delimiter='\t') if d['gene_id']==partner)
 print('CURRENT',json.dumps({k:d[k] for k in ['discovery_within_r','unadjusted_within_r','atlas_within_r','validation_directional_p','validation_directional_BH','replicated','partial_r_beyond_signature','atlas_partial_r']}))
# Examine the local gene-end evidence without assuming all downstream reads belong to either gene.
ends=[]
for p in sorted((root/'inputs/qc/three_prime').glob('*.gene_end_reads.tsv')):
 dd=rows(p)
 if len(dd)!=2:continue
 z={'sample_uid':p.name.split('.')[0]}
 for label,g in [('focal',focal),('partner',partner)]:
  a=int(dd[g]['inside']);b=int(dd[g]['downstream_1_to_W']);z[label+'_inside']=a;z[label+'_downstream']=b;z[label+'_downstream_fraction']=b/(a+b) if a+b else None
 ends.append(z)
if ends:
 with open(out/'gene_end_evidence.tsv','w') as f:
  w=csv.DictWriter(f,fieldnames=list(ends[0]),delimiter='\t');w.writeheader();w.writerows(ends)
 print('END_QC',len(ends),{k:float(np.median([d[k] for d in ends])) for k in ['focal_downstream_fraction','partner_downstream_fraction']})
