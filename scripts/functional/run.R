#!/usr/bin/env Rscript
# Retrospective, homology-transferred functional analysis of a fixed co-expression result.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(ggplot2);library(patchwork)})
setDTthreads(2);set.seed(20261109)
ROOT<-normalizePath('.');O<-file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "functional");A<-"inputs/functional_annotations";dir.create(O,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path('results','.cache'),recursive=TRUE,showWarnings=FALSE)
Sys.setenv(XDG_CACHE_HOME=file.path(ROOT,'results','.cache'))
NPERM<-10000L;FOCAL<-'Vitvi05_01chr01g05750'
readtab<-function(f,...)if(grepl('gz$',f))fread(cmd=paste('gzip -dc',shQuote(f)),...)else fread(f,...)
save<-function(d,f)fwrite(d,file.path(O,f),sep='\t',na='NA')
msg<-function(...){cat(format(Sys.time(),'%H:%M:%S'),...,'\n');flush.console()}
x<-readRDS(file.path(file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression"),"analysis_matrices.rds"))
d<-readtab(file.path(file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression"),"discovery_all_genes.tsv.gz"))
v<-readtab(file.path(file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression"),"atlas_all_genes.tsv.gz"))
rep<-readtab(file.path(file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression"),"replicated_partners.tsv"))
sig<-readtab('inputs/embryogenic_signature/signature_genes.tsv')[class=='up',gene_id]
mapping<-readtab(file.path(A,'gene_to_agi.tsv'));stopifnot(!anyDuplicated(mapping$gene_id))
go<-readtab(file.path(A,'agi_to_go.tsv.gz'));kegg<-readtab(file.path(A,'agi_to_kegg.tsv'))
terms<-rbind(readtab(file.path(A,'go_terms.tsv')),readtab(file.path(A,'kegg_terms.tsv')))
atann<-unique(rbind(go,kegg));ann<-unique(merge(mapping,atann,by='at_gene',allow.cartesian=TRUE))
common<-setdiff(intersect(x$discovery_genes,x$validation_genes),FOCAL)
pos<-d[stable_discovery==TRUE & discovery_within_r>0 & gene_id %in% common,gene_id]
repgenes<-rep$gene_id;outside<-setdiff(repgenes,sig)
stopifnot(length(repgenes)==106,length(pos)==521,length(outside)==54,all(repgenes %in% pos))
ann<-ann[gene_id %in% common]
save(ann,'transferred_gene_annotations.tsv.gz')
sets<-list(all_common=common,positive_candidates=pos,replicated=repgenes,replicated_outside_signature=outside)
coverage<-rbindlist(lapply(names(sets),function(n)rbindlist(lapply(c('BP','MF','CC','KEGG'),function(o){
 ids<-sets[[n]];g<-ann[ontology==o & gene_id %in% ids]
 data.table(set=n,ontology=o,n_genes=length(ids),n_with_homolog=sum(ids %in% mapping$gene_id),
 n_annotated=uniqueN(g$gene_id),n_unique_homologs=uniqueN(g$at_gene))}))))
save(coverage,'annotation_coverage.tsv');print(coverage)

# Fixed term family: sizes 10..500 in the shared, annotated gene universe.
# Same family is tested for both backgrounds and in collapsed-homolog sensitivity.
sizes<-ann[,.(K=uniqueN(gene_id)),by=.(ontology,term_id)]
eligible<-sizes[K>=10 & K<=500];save(merge(eligible,terms,by=c('term_id','ontology')),'eligible_terms.tsv')
comparisons<-list(
 replicated_vs_common=list(set=repgenes,background=common),
 replicated_vs_positive_candidates=list(set=repgenes,background=pos),
 outside_signature_vs_common=list(set=outside,background=setdiff(common,sig)),
 outside_signature_vs_positive_candidates=list(set=outside,background=setdiff(pos,sig)))

ora<-function(ids,background,o,unit='gene'){
 an<-ann[ontology==o];allanno<-unique(an$gene_id)
 ids<-intersect(ids,background);bg<-intersect(background,allanno);sel<-intersect(ids,allanno)
 if(unit=='homolog'){
  bg<-unique(mapping[gene_id %in% bg,at_gene]);sel<-unique(mapping[gene_id %in% sel,at_gene]);
  aa<-unique(an[,.(item=at_gene,term_id)])
 }else aa<-unique(an[,.(item=gene_id,term_id)])
 ee<-eligible[ontology==o,.(term_id)];N<-length(bg);n<-length(sel)
 bgsets<-split(aa[item %in% bg]$item,aa[item %in% bg]$term_id)
 selsets<-split(aa[item %in% sel]$item,aa[item %in% sel]$term_id)
 ee[,K:=vapply(term_id,function(t)length(unique(bgsets[[t]])),integer(1))]
 ee[,k:=vapply(term_id,function(t)length(unique(selsets[[t]])),integer(1))]
 ee[,`:=`(ontology=o,unit=unit,N_annotated=N,n_annotated=n,N_before_annotation=length(background),n_before_annotation=length(ids))]
 ee[,expected:=n*K/max(N,1)]
 ee[,fold_enrichment:=fifelse(expected>0,k/expected,NA_real_)]
 ee[,p:=if(N>0 && n>0)phyper(k-1,K,N-K,n,lower.tail=FALSE) else 1]
 ee[,BH:=p.adjust(p,'BH')]
 aaGenes<-an[gene_id %in% intersect(ids,background)]
 ee[,supporting_grape_genes:=vapply(term_id,function(t)paste(sort(unique(aaGenes[term_id==t,gene_id])),collapse=';'),character(1))]
 ee[,supporting_homologs:=vapply(term_id,function(t)paste(sort(unique(aaGenes[term_id==t,at_gene])),collapse=';'),character(1))]
 merge(ee,terms,by=c('term_id','ontology'),all.x=TRUE)
}
msg('Over-representation analysis')
oras<-rbindlist(lapply(names(comparisons),function(n){c<-comparisons[[n]]
 rbindlist(lapply(c('BP','MF','CC','KEGG'),function(o)rbindlist(lapply(c('gene','homolog'),function(u){
 a<-ora(c$set,c$background,o,u);a[,comparison:=n];a}))))}))
setorder(oras,comparison,ontology,unit,BH,p)
save(oras,'overrepresentation_all_tests.tsv.gz');save(oras[BH<.05],'overrepresentation_significant.tsv')
msg('Significant ORA terms:');print(oras[BH<.05,.(terms=.N),by=.(comparison,ontology,unit)])

# Rank-based competitive tests use all annotated genes shared between experiments.
# Statistic: probability a term gene ranks above a non-term gene by focal correlation (AUC).
# Permuting the focal vector within groups keeps the entire gene-gene covariance matrix intact.
unitrows<-function(M)M/pmax(sqrt(rowSums(M*M)),1e-12)
rank_test<-function(o,context){
 an<-unique(ann[ontology==o & term_id %in% eligible[ontology==o,term_id],.(gene_id,term_id)])
 # Background includes every gene with any annotation in this ontology, including those
 # whose own annotations are outside the selected term-size range.
 ids<-sort(unique(ann[ontology==o,gene_id]));N<-length(ids)
 ee<-eligible[ontology==o,term_id];an<-an[gene_id %in% ids]
 membership<-sparseMatrix(i=match(an$gene_id,ids),j=match(an$term_id,ee),x=1,dims=c(N,length(ee)))
 m<-as.numeric(colSums(membership));stopifnot(all(m>=10 & m<=500),all(m<N))
 if(context=='discovery'){M<-x$discovery_within;genes<-x$discovery_genes;g<-x$discovery_groups}
 else{M<-x$validation_within;genes<-x$validation_genes;g<-x$validation_groups}
 focal<-M[match(FOCAL,genes),];R<-unitrows(M[match(ids,genes),,drop=FALSE]);focal<-focal/sqrt(sum(focal*focal))
 obsr<-as.numeric(R%*%focal);obsrank<-rank(obsr,ties.method='average')
 rankmean<-as.numeric(crossprod(membership,obsrank))/m
 obs<-(rankmean-(m+1)/2)/(N-m)
 meanr<-as.numeric(crossprod(membership,obsr))/m
 blocks<-split(seq_along(g),g);tail<-integer(length(ee));nulsum<-numeric(length(ee));nulsq<-nulsum
 msg('Rank test',o,context,'genes',N,'terms',length(ee),'permutations',NPERM)
 for(st in seq(1,NPERM,by=100)){
  b<-min(100,NPERM-st+1)
  P<-replicate(b,{f<-focal;for(j in blocks)f[j]<-sample(focal[j]);f})
  rr<-R%*%P;rankmat<-apply(rr,2,rank,ties.method='average')
  means<-as.matrix(crossprod(membership,rankmat))/m
  auc<-(means-(m+1)/2)/(N-m)
  tail<-tail+rowSums(auc>=obs-1e-12);nulsum<-nulsum+rowSums(auc);nulsq<-nulsq+rowSums(auc*auc)
 }
 p<-(tail+1)/(NPERM+1);save(data.table(gene_id=ids,r=obsr,rank=obsrank),paste0('ranked_genes_',o,'_',context,'.tsv.gz'))
 data.table(term_id=ee,ontology=o,context=context,N=N,K=m,AUC=obs,mean_r=meanr,
 p=p,BH=p.adjust(p,'BH'),null_mean=nulsum/NPERM,null_sd=sqrt(pmax(0,nulsq/NPERM-(nulsum/NPERM)^2)))
}
ranked<-rbindlist(lapply(c('BP','KEGG'),function(o)rbind(rank_test(o,'discovery'),rank_test(o,'validation'))))
save(merge(ranked,terms,by=c('term_id','ontology')),'ranked_gene_set_all_tests.tsv')
dd<-ranked[context=='discovery'];vv<-ranked[context=='validation']
joint<-merge(dd,vv,by=c('term_id','ontology'),suffixes=c('_discovery','_validation'))
joint[,p_conjunction:=pmax(p_discovery,p_validation)]
joint[,BH_conjunction:=p.adjust(p_conjunction,'BH'),by=ontology]
joint[,reproduced_positive:=BH_conjunction<.05 & AUC_discovery>.5 & AUC_validation>.5 & mean_r_discovery>0 & mean_r_validation>0]
joint<-merge(joint,terms,by=c('term_id','ontology'))
joint[,supporting_replicated_genes:=vapply(term_id,function(t)paste(sort(unique(ann[term_id==t & gene_id %in% repgenes,gene_id])),collapse=';'),character(1))]
setorder(joint,ontology,BH_conjunction,p_conjunction)
save(joint,'ranked_gene_sets_conjunction.tsv');save(joint[reproduced_positive==TRUE],'ranked_gene_sets_reproduced.tsv')
msg('Ranked processes supported in both datasets:');print(joint[,.(tested=.N,reproduced=sum(reproduced_positive)),by=ontology])

# Summary of the four ORA contrasts and the homolog-collapse sensitivity.
ora_join<-merge(oras[unit=='gene'],oras[unit=='homolog',.(comparison,ontology,term_id,homolog_k=k,homolog_BH=BH,homolog_fold=fold_enrichment)],by=c('comparison','ontology','term_id'))
ora_join[,survives_homolog_collapse:=BH<.05 & homolog_BH<.05]
save(ora_join,'enrichment_with_homolog_sensitivity.tsv.gz')
save(ora_join[BH<.05],'enrichment_significant_with_homolog_sensitivity.tsv')

# Greedy reduction for display only: discard terms with >=70% Jaccard overlap
# with an already displayed term. All tested terms remain in the full tables.
reduce_terms<-function(ids,maxn=10){
 selected<-character();members<-split(ann$gene_id,ann$term_id);members<-lapply(members,unique)
 for(t in ids){j<-if(length(selected))vapply(selected,function(s)length(intersect(members[[s]],members[[t]]))/length(union(members[[s]],members[[t]])),numeric(1)) else 0
  if(max(j)<.70)selected<-c(selected,t)
  if(length(selected)>=maxn)break
 };selected
}
top<-oras[comparison=='replicated_vs_common' & ontology=='BP' & unit=='gene' & BH<.05]
top<-top[term_id %in% reduce_terms(top$term_id,10)]
disp<-ora_join[comparison %in% c('replicated_vs_common','replicated_vs_positive_candidates') & ontology=='BP' & unit=='gene' & term_id %in% top$term_id]
disp[,label:=vapply(term_name,function(s)paste(strwrap(s,43),collapse='\n'),character(1))]
disp[,comparison_label:=fifelse(comparison=='replicated_vs_common','Shared measurable genes','Positive discovery candidates')]
disp[,label:=factor(label,levels=rev(unique(top[,vapply(term_name,function(s)paste(strwrap(s,43),collapse='\n'),character(1))])))]
theme_set(theme_classic(base_size=10)+theme(legend.position='bottom',plot.title=element_text(face='bold')))
p1<-ggplot(disp,aes(fold_enrichment,label))+geom_point(aes(size=k,colour=-log10(pmax(BH,1e-300)),shape=BH<.05))+
 scale_colour_viridis_c(option='C')+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16))+
 facet_wrap(~comparison_label)+labs(title='A  GO biological processes of replicated partners',x='Fold enrichment',y=NULL,size='Grapevine genes',colour='-log10(BH)',shape='BH < 0.05')
jtop<-joint[ontology=='BP' & reproduced_positive==TRUE];jtop<-jtop[term_id %in% reduce_terms(jtop$term_id,10)]
rd<-rbind(jtop[,.(term_id,term_name,dataset='Discovery',AUC=AUC_discovery)],jtop[,.(term_id,term_name,dataset='Independent callus',AUC=AUC_validation)])
rd[,label:=vapply(term_name,function(s)paste(strwrap(s,43),collapse='\n'),character(1))]
rd[,label:=factor(label,levels=rev(unique(jtop[,vapply(term_name,function(s)paste(strwrap(s,43),collapse='\n'),character(1))])))]
p2<-ggplot(rd,aes(AUC,label,colour=dataset))+geom_vline(xintercept=.5,linetype=2,colour='grey60')+geom_point(size=2.5,position=position_dodge(width=.4))+
 scale_colour_manual(values=c('Discovery'='#2678ae','Independent callus'='#c45236'))+
 labs(title='B  Rank-based processes supported in both datasets',x='Probability a process gene ranks above a background gene',y=NULL,colour=NULL)
if(!nrow(disp))p1<-ggplot()+annotate('text',x=0,y=0,label='No GO BP over-representation at BH < 0.05')+theme_void()
if(!nrow(rd))p2<-ggplot()+annotate('text',x=0,y=0,label='No GO BP rank-based conjunction at BH < 0.05')+theme_void()
fig<-p1/p2+plot_layout(heights=c(1.2,1))+plot_annotation(caption='Exploratory annotation transfer through existing Arabidopsis best-hit homologs. GO terms are reduced for display only.\nRank tests use within-group permutations; conjunction requires support in both datasets. Complete results and homolog-collapse sensitivity are supplied.')
for(ext in c('pdf','png'))ggsave(file.path(O,paste0('functional_summary.',ext)),fig,width=14,height=12,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else 'png')
save(disp,'displayed_GO_enrichment.tsv');save(jtop,'displayed_ranked_processes.tsv')
capture.output(sessionInfo(),file=file.path(O,'sessionInfo.txt'))
saveRDS(list(common=common,positive_candidates=pos,replicated=repgenes,outside_signature=outside,
 mapping=mapping,annotations=ann,eligible_terms=eligible,ora=oras,ranked=ranked,joint=joint,
 parameters=list(seed=20261109,NPERM=NPERM,min_term_size=10,max_term_size=500)),file.path(O,'functional_analysis.rds'))
source('scripts/functional/plot.R')
msg('DONE',file.path(ROOT,O))
