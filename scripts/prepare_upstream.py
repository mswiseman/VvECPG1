#!/usr/bin/env python3
"""Stage a writable original-layout project. Frozen package inputs are never edited."""
from pathlib import Path
import csv,gzip,shutil,os
P=Path(__file__).resolve().parents[1]; O=P/'results/upstream_workspace'
if (O/'.prepared').exists():
 print('Already staged:',O,'(preserving existing results)');raise SystemExit(0)
O.mkdir(parents=True,exist_ok=True)
A=P/'source_addenda/A_upstream_workflow';B=P/'source_addenda/B_atlas_salmon';C=P/'source_addenda/C_de_tables'
def cp(src,dst):
 dst=O/dst;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
def tree(src,dst):shutil.copytree(src,O/dst,dirs_exist_ok=True)
def unzip(src,dst):
 dst=O/dst;dst.parent.mkdir(parents=True,exist_ok=True)
 with gzip.open(src,'rb') as r,dst.open('wb') as w:shutil.copyfileobj(r,w)
tree(P/'upstream/scripts',Path('scripts'));tree(P/'upstream/config',Path('config'));tree(P/'upstream/workflow',Path('workflow'));tree(P/'upstream/envs',Path('envs'))
cp(P/'upstream/Makefile',Path('Makefile'))
cp(P/'inputs/metadata/source_worksheet.xlsx',Path('Samples information for RNA-seq data.xlsx'))
for f in ('samples.tsv','sample_manifest.tsv','fastq_to_sample_uid.tsv'):cp(A/'checksums'/f,Path('metadata')/f)
cp(P/'inputs/counts/pn_v51.primary.counts.tsv',Path('results/pn_v51/counts/pn_v51.primary.counts.tsv'))
for f in (P/'inputs/counts').glob('*'):cp(f,Path('results/pn40024_v5.1/counts')/f.name)
cp(P/'inputs/pn_v51.per_sample_qc.tsv',Path('results/pn_v51/qc/pn_v51.per_sample_qc.tsv'))
for f in (P/'inputs/focal_network').glob('*'):cp(f,Path('results/pn40024_v5.1/focal_network')/f.name)
for f in (P/'inputs/embryogenic_signature').glob('*'):cp(f,Path('results/sra_atlas/embryogenic_signature')/f.name)
for f in (P/'inputs/atlas').glob('*'):cp(f,Path('results/sra_atlas')/f.name)
for f in C.glob('*.tsv'):cp(f,Path('results/pn40024_v5.1/de')/f.name)
for f in (C/'deseq2_objects').glob('*'):cp(f,Path('results/pn40024_v5.1/de')/f.name)
for f in (C/'de_tables').glob('*.gz'):unzip(f,Path('results/pn40024_v5.1/de/tables')/f.name.removesuffix('.gz'))
for f in (C/'focal').glob('*'):cp(f,Path('results/pn40024_v5.1/target_gene')/f.name)
unzip(B/'salmon_index/tx2gene.tsv.gz',Path('references/salmon_v51/tx2gene.tsv'))
# Runtime relative symlink avoids making a second 480-MB quantification copy. Read-only by convention.
q=O/'results/sra_atlas/quant'
if not q.exists():q.symlink_to(os.path.relpath(B/'quant',q.parent),target_is_directory=True)
# Recover annotation-only lookup columns from an archived all-gene DE table.
# This affects labels, not the expression matrix or statistical calculations.
f=next((C/'de_tables').glob('*.gz'))
cols=['symbol','descriptor','vendor_ensembl_r62_ids','vendor_mapping_status','vitvi_v4_ids','vitvi_v3_ids','VIT_v1_ids','max_utr3_bp']
dst=O/'results/crosswalk/gene_id_lookup_wide.tsv';dst.parent.mkdir(parents=True,exist_ok=True)
with gzip.open(f,'rt') as h,dst.open('w') as out:
 r=csv.DictReader(h,delimiter='\t');w=csv.DictWriter(out,fieldnames=['pn_v5_1_gene_id']+cols,delimiter='\t');w.writeheader()
 for x in r:w.writerow(dict(pn_v5_1_gene_id=x['gene_id'],**{c:x.get(c,'') for c in cols}))
# Fresh atlas quantification uses the final curated cohort, including separately flagged sensitivity runs.
with (P/'inputs/atlas/focal_atlas_per_run.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
sel=[x for x in rows if x['qc_ok'].lower()=='true' and x['in_selection'].lower()=='true' and x['tissue_class']]
md=O/'results/sra_atlas/metadata';md.mkdir(parents=True,exist_ok=True)
with (md/'selection.tsv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(sel)
(md/'run_list.txt').write_text(''.join(x['run_accession']+'\n' for x in sel))
for name in ('curated_studies.tsv','run_overrides.tsv'):cp(P/'inputs/atlas'/name,Path('results/sra_atlas/metadata')/name)
# Frozen homology outputs are copied to writable locations before any rescan.
E=P/'source_addenda/E_vvecpg1_homology'
for f in E.iterdir():
 if f.name in ('scripts','SHA256SUMS','README.md','ANALYSIS_NOTES.md'):continue
 if f.is_dir():tree(f,Path('results/homology')/f.name)
 else:cp(f,Path('results/homology')/f.name)
(O/'.prepared').write_text('Staged from packaged immutable inputs; all generated results stay here.\n')
print('Staged writable upstream project:',O)
print('Raw FASTQ/BAMs, reference genomes/index binaries and tool environments remain external; see upstream/README.md.')
