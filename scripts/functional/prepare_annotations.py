"""Prepare auditable homology-transferred GO/KEGG annotations; Python standard library only."""
from pathlib import Path
import csv, gzip, re, hashlib, json
from collections import defaultdict, Counter
from functools import lru_cache

ROOT=Path(__file__).resolve().parents[2]
A=ROOT/'inputs/functional_annotations'
sources={
 'tair.gaf.gz':'https://current.geneontology.org/annotations/tair.gaf.gz',
 'go-basic.obo':'https://current.geneontology.org/ontology/go-basic.obo',
 'kegg_ath_gene_pathway.tsv':'https://rest.kegg.jp/link/pathway/ath',
 'kegg_ath_pathway_names.tsv':'https://rest.kegg.jp/list/pathway/ath'}
def write(path,fields,rows):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'wt',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(fields);w.writerows(rows)

with gzip.open(ROOT/'inputs/focal_network/focal_coexpression_all_genes.tsv.gz','rt') as f:
 mapping={r['gene_id']:r['at_gene'].upper().split('.')[0] for r in csv.DictReader(f,delimiter='\t')
          if re.fullmatch(r'AT[1-5MC]G[0-9]{5}',r['at_gene'].upper().split('.')[0])}
write(A/'gene_to_agi.tsv',['gene_id','at_gene'],sorted(mapping.items()))
wanted=set(mapping.values())
terms={};aliases={};cur=None;version=''
def store(t):
 if not t or 'id' not in t:return
 if not t.get('obsolete',False):
  terms[t['id']]=t
  for a in t['alt']:aliases[a]=t['id']
with (A/'go-basic.obo').open() as f:
 for raw in f:
  line=raw.strip()
  if line.startswith('data-version:'):version=line.split(': ',1)[1]
  if line=='[Term]':store(cur);cur={'parents':[],'alt':[]}
  elif line.startswith('['):store(cur);cur=None
  elif cur is not None:
   if line.startswith('id: '):cur['id']=line[4:]
   elif line.startswith('name: '):cur['name']=line[6:]
   elif line.startswith('namespace: '):cur['namespace']=line[11:]
   elif line.startswith('is_a: '):cur['parents'].append(line.split()[1])
   elif line.startswith('relationship: part_of '):cur['parents'].append(line.split()[2])
   elif line.startswith('alt_id: '):cur['alt'].append(line.split()[1])
   elif line=='is_obsolete: true':cur['obsolete']=True
store(cur)
@lru_cache(None)
def ancestors(g):
 t=terms[g];out={g}
 for p in t['parents']:
  p=aliases.get(p,p)
  if p in terms and terms[p]['namespace']==t['namespace']:out.update(ancestors(p))
 return frozenset(out)

annotations=defaultdict(set);evidence=Counter();stats=Counter();headers=[];current_symbols={}
with gzip.open(A/'tair.gaf.gz','rt') as f:
 for line in f:
  if line.startswith('!'):
   if len(headers)<30:headers.append(line.rstrip())
   continue
  r=line.rstrip('\n').split('\t')
  if len(r)<15:continue
  ag=r[1].upper().split('.')[0]
  if ag not in wanted:continue
  current_symbols.setdefault(ag,(r[2],r[9]))
  stats['mapped_homolog_direct_rows']+=1
  if 'NOT' in r[3].split('|'):stats['excluded_NOT']+=1;continue
  if r[6]=='ND':stats['excluded_ND']+=1;continue
  go=aliases.get(r[4],r[4])
  if go not in terms:stats['obsolete_or_missing_GO']+=1;continue
  annotations[ag].update(ancestors(go));evidence[r[6]]+=1
ns={'biological_process':'BP','molecular_function':'MF','cellular_component':'CC'}
write(A/'agi_to_go.tsv.gz',['at_gene','term_id','ontology'],
 ((ag,g,ns[terms[g]['namespace']]) for ag,gs in sorted(annotations.items()) for g in sorted(gs)))
write(A/'go_terms.tsv',['term_id','term_name','ontology'],
 ((g,t['name'],ns[t['namespace']]) for g,t in sorted(terms.items()) if t['namespace'] in ns))
write(A/'go_parent_edges.tsv',['child','parent'],
 ((g,aliases.get(p,p)) for g,t in sorted(terms.items()) for p in t['parents'] if aliases.get(p,p) in terms))
write(A/'tair_current_gene_names.tsv',['at_gene','tair_symbol','tair_description'],
 ((ag,*nm) for ag,nm in sorted(current_symbols.items())))

names={}
with (A/'kegg_ath_pathway_names.tsv').open() as f:
 for row in csv.reader(f,delimiter='\t'):
  if len(row)>=2:names[row[0].replace('path:','')]=row[1].removesuffix(' - Arabidopsis thaliana (thale cress)')
kegg=set()
with (A/'kegg_ath_gene_pathway.tsv').open() as f:
 for row in csv.reader(f,delimiter='\t'):
  if len(row)!=2:continue
  ag=row[0].replace('ath:','').upper();p=row[1].replace('path:','')
  if ag in wanted and not p.startswith(('ath011','ath012')):kegg.add((ag,p,'KEGG'))
write(A/'agi_to_kegg.tsv',['at_gene','term_id','ontology'],sorted(kegg))
write(A/'kegg_terms.tsv',['term_id','term_name','ontology'],
 ((p,n,'KEGG') for p,n in sorted(names.items()) if not p.startswith(('ath011','ath012'))))
metadata={'download_date':'2026-11-08','GO_ontology_version':version,'GAF_header':headers,
 'source_mapping':'Existing focal_coexpression_all_genes.tsv.gz at_gene best-hit homology labels',
 'mapping_limitation':'Underlying DIAMOND e-values, coverages and RBH flags were not supplied; cannot reapply stricter upstream thresholds.',
 'GO_propagation':'is_a and part_of only, within ontology; NOT and ND excluded; nonobsolete terms, alt IDs resolved',
 'KEGG':'Arabidopsis pathways transferred through provided best-hit homology; NOT direct Vitis pathway annotation. Overview maps ath011* and ath012* excluded.',
 'n_mapped_grape_genes':len(mapping),'n_unique_homologs':len(wanted),
 'n_homologs_GO':len(annotations),'GAF_filter_counts':dict(stats),'GAF_evidence_counts':dict(evidence)}
(OUT/'annotation_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
write(OUT/'annotation_sources.tsv',['filename','url','download_date','bytes','sha256'],
 ((name,url,'2026-11-08',(A/name).stat().st_size,hashlib.sha256((A/name).read_bytes()).hexdigest()) for name,url in sources.items()))
print(json.dumps({k:metadata[k] for k in ['GO_ontology_version','n_mapped_grape_genes','n_unique_homologs','n_homologs_GO']},indent=2))
