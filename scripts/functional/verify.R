#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(Matrix)})
O<-file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "functional")
x<-readRDS(file.path(O,'functional_analysis.rds'));checks<-list()
check<-function(n,v){stopifnot(isTRUE(v));checks[[length(checks)+1]]<<-data.table(check=n,passed=TRUE)}
a<-x$ora;an<-x$annotations
check('Fixed replicated gene set has 106 genes',length(x$replicated)==106)
check('Positive-candidate background has 521 genes',length(x$positive_candidates)==521)
check('Non-signature subset has 54 genes',length(x$outside_signature)==54)
check('Gene annotations unique by homolog/gene/term',!anyDuplicated(an[,.(gene_id,at_gene,term_id)]))
check('Every annotation is in the shared measured universe',all(an$gene_id %in% x$common))
check('Every target gene belongs to positive background',all(x$replicated %in% x$positive_candidates))
check('Valid contingency table counts',all(a$k<=a$n_annotated & a$k<=a$K & a$n_annotated<=a$N_annotated & a$K<=a$N_annotated))
families<-split(seq_len(nrow(a)),paste(a$comparison,a$ontology,a$unit))
check('All 32 planned enrichment families present',length(families)==32)
for(n in names(families)){
 q<-a[families[[n]]]
 check(paste('BH includes all terms',n),max(abs(p.adjust(q$p,'BH')-q$BH))<1e-12)
}
sel<-a[unit=='gene' & comparison=='replicated_vs_common' & ontology=='BP' & BH<.05]
for(i in seq_len(nrow(sel))){
 t<-sel[i];z<-matrix(c(t$k,t$n_annotated-t$k,t$K-t$k,t$N_annotated-t$K-t$n_annotated+t$k),2,byrow=TRUE)
 check(paste('Hypergeometric agrees with exact Fisher',t$term_id),abs(fisher.test(z,alternative='greater')$p.value-t$p)<1e-12)
 check(paste('Supporting gene count',t$term_id),length(strsplit(t$supporting_grape_genes,';',fixed=TRUE)[[1]])==t$k)
}
j<-x$joint
check('Conjunction p is maximum of dataset p values',max(abs(j$p_conjunction-pmax(j$p_discovery,j$p_validation)))<1e-12)
for(o in unique(j$ontology))check(paste('Conjunction BH',o),max(abs(p.adjust(j[ontology==o,p_conjunction],'BH')-j[ontology==o,BH_conjunction]))<1e-12)
for(o in c('BP','KEGG'))for(ctx in c('discovery','validation')){
 r<-fread(cmd=paste('gzip -dc',shQuote(file.path(O,paste0('ranked_genes_',o,'_',ctx,'.tsv.gz')))))
 test<-x$ranked[ontology==o & context==ctx][1]
 ids<-an[ontology==o & term_id==test$term_id,gene_id];q<-r[gene_id %in% ids]
 auc<-(mean(q$rank)-(nrow(q)+1)/2)/(nrow(r)-nrow(q))
 check(paste('Observed rank AUC independently reconstructed',o,ctx),abs(auc-test$AUC)<1e-12)
 check(paste('Observed set mean r',o,ctx),abs(mean(q$r)-test$mean_r)<1e-12)
}
checks<-rbindlist(checks);dir.create('results/validation',recursive=TRUE,showWarnings=FALSE);fwrite(checks,'results/validation/functional_checks.tsv',sep='\t');cat(nrow(checks),'functional-analysis checks passed\n')
