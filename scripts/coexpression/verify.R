#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
p<-paste0(file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression"), "/")
x<-readRDS(paste0(p,'analysis_matrices.rds'))
tab<-fread(paste0(p,'fixed_candidates_validation.tsv'))
rp<-tab[replicated==TRUE]
sig<-fread('inputs/embryogenic_signature/signature_genes.tsv')[class=='up',gene_id]
focal<-'Vitvi05_01chr01g05750'
W<-x$discovery_within;genes<-x$discovery_genes;fi<-match(focal,genes)
checks<-list()
check<-function(name,value){stopifnot(isTRUE(value));checks[[length(checks)+1]]<<-data.table(check=name,passed=TRUE)}
check('90 treated discovery samples / 30 groups',ncol(W)==90 && uniqueN(x$discovery_groups)==30)
check('27 independent validation samples / 9 groups',ncol(x$validation_within)==27 && uniqueN(x$validation_groups)==9)
check('Focal excluded from signature',!focal %in% sig)
check('Within-group matrix rows have zero group means',max(abs(sapply(unique(x$discovery_groups),function(g)rowMeans(W[,x$discovery_groups==g,drop=FALSE]))))<1e-10)
for(k in head(which(genes %in% tab$gene_id),12)){
 id<-genes[k];t<-tab[gene_id==id]
 check(paste('Direct Pearson r',id),abs(cor(W[fi,],W[k,])-t$discovery_within_r)<1e-10)
 up<-setdiff(intersect(sig,genes),id);score<-colMeans(W[match(up,genes),,drop=FALSE])
 pr<-cor(resid(lm(W[fi,]~score)),resid(lm(W[k,]~score)))
 check(paste('Leave-partner-out partial r',id),abs(pr-t$partial_r_beyond_signature)<1e-10)
 loo<-sapply(seq_len(ncol(W)),function(i){m<-W[,-i,drop=FALSE];g<-x$discovery_groups[-i];for(a in unique(g))m[,g==a]<-m[,g==a,drop=FALSE]-rowMeans(m[,g==a,drop=FALSE]);cor(m[fi,],m[k,])})
 check(paste('Leave-one-sample recomputes group means',id),abs(min(loo*sign(t$discovery_within_r))-t$minimum_signed_leave_one_sample_out_r)<1e-10)
}
ok<-!is.na(tab$validation_directional_p)
check('Validation BH is over measured fixed candidates',max(abs(p.adjust(tab$validation_directional_p[ok],'BH')-tab$validation_directional_BH[ok]))<1e-10)
check('Every replicated partner has matching sign and BH < 0.05',all(rp$same_sign & rp$validation_directional_BH<.05))
check('Replicated count',nrow(rp)==106)
check('Independent sample identities unique',!anyDuplicated(fread(paste0(p,'validation_samples.tsv'))$sample_accession))
check('Signature all-library r reproduced',abs(fread(paste0(p,'signature_reproduction.tsv'))[context=='all_libraries' & transformation=='adjusted',within_r]-.7382236)<1e-6)
check('Legacy overlap matches the report',fread(paste0(p,'historical_overlap_check.tsv'))[method=='legacy_fixed_cutoff_check',overlap]==11)
checks<-rbindlist(checks);dir.create('results/validation',recursive=TRUE,showWarnings=FALSE);fwrite(checks,'results/validation/coexpression_checks.tsv',sep='\t');cat(nrow(checks),'scientific/numerical checks passed\n')
