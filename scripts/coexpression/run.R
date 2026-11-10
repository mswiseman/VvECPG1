#!/usr/bin/env Rscript
# Local, reproducible focal-gene follow-up. No upstream files are modified.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggrepel)})
setDTthreads(2)
set.seed(20261108)
ROOT <- normalizePath('.')
dir.create(file.path(ROOT,'results','.cache'),recursive=TRUE,showWarnings=FALSE)
Sys.setenv(XDG_CACHE_HOME=file.path(ROOT,'results','.cache'))
OUT <- file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
FOCAL <- 'Vitvi05_01chr01g05750'
NPERM <- 1000L; NVALID <- 10000L; NBOOT <- 200L
logmsg <- function(...) {cat(format(Sys.time(),'%H:%M:%S'),...,'\n');flush.console()}
readtab <- function(f,...) if(grepl('\\.gz$',f)) fread(cmd=paste('gzip -dc',shQuote(f)),...) else fread(f,...)
save <- function(x,n) fwrite(x,file.path(OUT,n),sep='\t',na='NA')
row_unit <- function(M) M/pmax(sqrt(rowSums(M^2)),1e-12)
center <- function(M,g) {for(k in unique(g)){j<-which(g==k);M[,j]<-M[,j,drop=FALSE]-rowMeans(M[,j,drop=FALSE])};M}
zrows <- function(M){a<-t(scale(t(M)));a[!is.finite(a)]<-0;a}
cr <- function(M,x) as.numeric(row_unit(M)%*%(x/sqrt(sum(x*x))))
gmean <- function(M,g) sapply(unique(g),function(k)rowMeans(M[,g==k,drop=FALSE]))
sig <- readtab('inputs/embryogenic_signature/signature_genes.tsv')
up <- sig[class=='up',gene_id];stopifnot(!FOCAL %in% up)
ann <- readtab('inputs/focal_network/focal_coexpression_all_genes.tsv.gz')
cc <- readRDS(Sys.getenv('VVECPG1_ADJUSTED_RDS','inputs/focal_network/context_C_adjusted.rds'));meta<-as.data.table(cc$s)
stopifnot(identical(colnames(cc$adj),meta$sample_uid),!anyDuplicated(meta$sample_uid))
blind<-readtab('inputs/counts/pn_v51.vst_blind_visualization_only.tsv')
BV<-as.matrix(blind[,-1]);rownames(BV)<-blind$gene_id;BV<-BV[,meta$sample_uid]
nc<-readtab('inputs/counts/pn_v51.normalized_counts_deseq2.tsv')
NC<-as.matrix(nc[,-1]);rownames(NC)<-nc$gene_id;NC<-NC[,meta$sample_uid]
counts<-readtab('inputs/counts/pn_v51.primary.counts.tsv')
CM<-as.matrix(counts[,meta$sample_uid,with=FALSE]);rownames(CM)<-counts$gene_id
stopifnot(all(CM>=0),all(CM==round(CM)),FOCAL %in% rownames(CM))
save(meta[,.(sample_uid,batch,group,hormone,dose_value,dose_unit,time_hours,is_control)],'experiment_samples.tsv')
logmsg('Inputs:',nrow(cc$adj),'transformed genes;',ncol(cc$adj),'libraries;',length(up),'positive signature genes')

# Reproduce the original signature result, then remove duplicate control preparations.
orig<-readtab('inputs/embryogenic_signature/sample_scores.tsv')
repro<-list()
for(n in c('all_libraries','treated_only','early_treated','24h_treated')){
 ix<-switch(n,all_libraries=rep(TRUE,nrow(orig)),treated_only=!orig$is_control,
            early_treated=!orig$is_control & orig$time_hours %in% c(6,12,24),
            `24h_treated`=!orig$is_control & orig$time_hours==24)
 for(v in c('adjusted','unadjusted')){
  a<-orig[ix];g<-paste(a$batch,a$group);x<-a[[paste0('focal_z_',v)]];y<-a[[paste0('score_up_',v)]]
  repro[[length(repro)+1]]<-data.table(context=n,transformation=v,n_libraries=nrow(a),n_groups=uniqueN(g),
    n_biological=if(n=='all_libraries')123L else nrow(a),within_r=cor(x-ave(x,g),y-ave(y,g)))
 }
}
repro<-rbindlist(repro);save(repro,'signature_reproduction.tsv');print(repro)
cs<-readtab('inputs/embryogenic_signature/contrast_scores.tsv')
contr<-rbindlist(lapply(c('all40','24h'),function(k){a<-if(k=='24h')cs[time_hours==24] else cs
 data.table(context=k,n_contrasts=nrow(a),pearson_raw_score=cor(a$focal_lfc,a$score_up),
 pearson_background_adjusted=cor(a$focal_lfc,a$score_up-a$score_up_null_median),
 spearman_background_adjusted=cor(a$focal_lfc,a$score_up-a$score_up_null_median,method='spearman'))}))
save(contr,'contrast_correlations.tsv')
save(cs[time_hours==24][,programme_shift:=score_up-score_up_null_median],'24h_contrast_scores.tsv')
loo24<-rbindlist(lapply(seq_len(nrow(cs[time_hours==24])),function(i){a<-cs[time_hours==24];b<-a[-i];data.table(left_out=a$contrast[i],r=cor(b$focal_lfc,b$score_up-b$score_up_null_median))}))
save(loo24,'24h_contrast_leave_one_out.tsv')

# Use the adjusted matrix exactly as supplied; filter expression before association testing.
early<-which(!meta$is_control & meta$time_hours %in% c(6,12,24) & meta$batch==1)
late<-which(!meta$is_control & meta$time_hours==168)
stopifnot(length(early)==90,length(late)==30,uniqueN(meta$group[early])==30)
keep<-rowSums(CM[,early,drop=FALSE]>=10)>=3
genes<-intersect(rownames(cc$adj),rownames(CM)[keep]);genes<-genes[startsWith(genes,'Vitvi05_01chr')]
stopifnot(FOCAL %in% genes)
E<-cc$adj[genes,early,drop=FALSE];g<-meta$group[early];Z<-zrows(E);W<-center(Z,g)
fi<-match(FOCAL,genes);x<-W[fi,];r<-cr(W,x)
U<-intersect(up,genes);ui<-match(U,genes);score<-colMeans(W[ui,,drop=FALSE])

# Pair-specific partial r: remove the tested partner from the positive-signature score.
partial_rows<-function(W,genes,up,fi){
 U<-intersect(up,genes);j<-match(U,genes);total<-colSums(W[j,,drop=FALSE]);n<-length(j)
 S<-matrix(total/n,nrow(W),ncol(W),byrow=TRUE)
 if(length(j))S[j,]<-(matrix(total,length(j),ncol(W),byrow=TRUE)-W[j,,drop=FALSE])/(n-1)
 ss<-rowSums(S*S);xx<-matrix(W[fi,],nrow(W),ncol(W),byrow=TRUE)
 xr<-xx-S*(rowSums(xx*S)/pmax(ss,1e-12));yr<-W-S*(rowSums(W*S)/pmax(ss,1e-12))
 list(r=rowSums(row_unit(xr)*row_unit(yr)),x=xr,y=yr,score=S)
}
part<-partial_rows(W,genes,up,fi)
Wu<-center(zrows(BV[genes,early,drop=FALSE]),g);ru<-cr(Wu,Wu[fi,])
G<-gmean(Z,g);G<-G-rowMeans(G);rb<-cr(G,G[fi,])
logmsg('Calling discovery partners using within-group permutations')

# Pooled empirical FDR, as in the original analysis. Separate positive/negative tails.
# This estimates a discovery-set error fraction, not per-gene adjusted p-values.
pool_fdr<-function(W,x,g,nperm=NPERM,grid=seq(.30,.95,.005)){
 R<-row_unit(W);obs<-cr(W,x);eligible<-seq_len(nrow(W));eligible<-eligible[rownames(W)!=FOCAL]
 obsP<-sapply(grid,function(t)sum(obs[eligible]>=t));obsN<-sapply(grid,function(t)sum(obs[eligible]<=-t))
 nullP<-numeric(length(grid));nullN<-nullP;tail<-integer(nrow(W))
 blocks<-split(seq_along(g),g)
 for(start in seq(1,nperm,by=100)){
  b<-min(100,nperm-start+1)
  P<-replicate(b,{xp<-x;for(j in blocks)xp[j]<-sample(x[j]);xp/sqrt(sum(xp*xp))})
  rr<-R%*%P;tail<-tail+rowSums(abs(rr)>=abs(obs)-1e-12)
  v<-sort(as.vector(rr[eligible,,drop=FALSE]));nv<-length(v)
  nullP<-nullP+(nv-findInterval(grid,v,left.open=TRUE));nullN<-nullN+findInterval(-grid,v)
 }
 fp<-nullP/nperm/pmax(obsP,1);fn<-nullN/nperm/pmax(obsN,1)
 tp<-if(any(fp<=.05 & obsP>0))grid[which(fp<=.05 & obsP>0)[1]] else Inf
 tn<-if(any(fn<=.05 & obsN>0))grid[which(fn<=.05 & obsN>0)[1]] else Inf
 list(r=obs,pos=tp,neg=tn,p=(tail+1)/(nperm+1),curve=data.table(abs_r=grid,observed_positive=obsP,
  observed_negative=obsN,expected_null_positive=nullP/nperm,expected_null_negative=nullN/nperm,
  pooled_fdr_positive=fp,pooled_fdr_negative=fn))
}
fd<-pool_fdr(W,x,g);save(fd$curve,'discovery_permutation_fdr.tsv')
called<-(r>=fd$pos | r<=-fd$neg);called[fi]<-FALSE
logmsg('Discovery thresholds +',fd$pos,'-',fd$neg,'called',sum(called))

# Stability: resample entire treatment groups (not individual residual values).
stab<-numeric(length(genes));signst<-numeric(length(genes));groups<-unique(g)
for(b in seq_len(NBOOT)){
 gs<-sample(groups,length(groups),replace=TRUE);ix<-unlist(lapply(gs,function(k)which(g==k)))
 br<-cr(W[,ix,drop=FALSE],x[ix]);stab<-stab+((br>=fd$pos & r>0)|(br<=-fd$neg & r<0))
 signst<-signst+(sign(br)==sign(r))
}
loo<-sapply(unique(meta$hormone[early]),function(h){ix<-which(meta$hormone[early]!=h);cr(W[,ix,drop=FALSE],x[ix])})
colnames(loo)<-paste0('r_without_',unique(meta$hormone[early]))
low<-!meta$sample_uid[early] %in% c('PPPJ8P_19','PPPJ8P_52','PPPJ8P_81')
# Recenter after removal: a two-replicate group otherwise retains a non-zero mean.
Wlow<-center(Z[,low,drop=FALSE],g[low]);rlow<-cr(Wlow,Wlow[fi,])
i24<-which(meta$time_hours[early]==24);r24<-cr(W[,i24,drop=FALSE],x[i24])
L<-zrows(cc$adj[genes,late,drop=FALSE]);WL<-center(L,paste(meta$batch[late],meta$group[late]));rlate<-cr(WL,WL[fi,])
stable<-called & stab/NBOOT>=.80 & apply(loo*sign(r),1,min)>=.30
# Exact leave-one-library sensitivity, re-centering the remaining two replicates.
# For a group of size n, removing one observation subtracts n/(n-1)*xi*yi
# from the centered sum of products (and the analogous sum of squares).
leave_one_min<-function(W,fi,g,sgn){
 xx<-W[fi,];xy<-as.numeric(W%*%xx);yy<-rowSums(W*W);sx<-sum(xx*xx)
 out<-rep(Inf,nrow(W));weights<-as.numeric(table(g)[g])/(as.numeric(table(g)[g])-1)
 for(i in seq_along(xx)){
  rr<-(xy-weights[i]*W[,i]*xx[i])/sqrt(pmax(yy-weights[i]*W[,i]^2,1e-12)*pmax(sx-weights[i]*xx[i]^2,1e-12))
  out<-pmin(out,rr*sgn)
 };out
}
leave_sample<-leave_one_min(W,fi,g,sign(r))
tab<-data.table(gene_id=genes,discovery_within_r=r,discovery_between_r=rb,unadjusted_within_r=ru,
  partial_r_beyond_signature=part$r,pooled_fdr_call=called,cluster_bootstrap_selection_frequency=stab/NBOOT,
  bootstrap_sign_frequency=signst/NBOOT,minimum_signed_leave_hormone_out_r=apply(loo*sign(r),1,min),
  r_without_flagged_libraries=rlow,r_24h=r24,r_7day=rlate,stable_discovery=stable,
  discovery_empirical_p=fd$p,in_positive_signature=genes %in% U,
  minimum_signed_leave_one_sample_out_r=leave_sample)
tab<-cbind(tab,as.data.table(loo));save(tab,'discovery_all_genes.tsv.gz')
save(meta[early,.(sample_uid,group,hormone,dose_value,dose_unit,time_hours,batch)],'discovery_samples.tsv')
save(data.table(n_samples=90,n_groups=30,n_genes=length(genes)-1,n_signature=length(U),
  positive_threshold=fd$pos,negative_threshold=fd$neg,n_called=sum(called),n_stable=sum(stable)),'discovery_summary.tsv')
logmsg('Stable discovery candidates:',sum(stable))

# Independent validation: only embryogenic callus, one row per biological BioSample.
runs<-readtab('inputs/atlas/focal_atlas_per_run.tsv')[study_accession=='PRJNA732451' & qc_ok==TRUE]
stopifnot(!anyDuplicated(runs$sample_accession),nrow(runs)==27)
runs[,condition:=sub('(_rep)?[0-9]+$','',sample_title)]
logmsg('Reading atlas expression matrix')
at<-readtab('inputs/atlas/gene_tpm.tsv.gz',select=c(1,match(runs$run_accession,names(readtab('inputs/atlas/gene_tpm.tsv.gz',nrows=0)))))
AG<-at[[1]];AT<-as.matrix(at[,-1]);rownames(AT)<-AG
eligible<-startsWith(AG,'Vitvi05_01chr') & rowMeans(AT>=5)>=.20
eligible[match(FOCAL,AG)]<-TRUE;AT<-AT[eligible,,drop=FALSE];ag<-rownames(AT);afi<-match(FOCAL,ag)
AZ<-zrows(log2(AT+1));AW<-center(AZ,runs$condition);ar<-cr(AW,AW[afi,]);ap<-partial_rows(AW,ag,up,afi)
save(runs,'validation_samples.tsv')
afd<-pool_fdr(AW,AW[afi,],runs$condition);save(afd$curve,'atlas_permutation_fdr.tsv')
av<-data.table(gene_id=ag,atlas_within_r=ar,atlas_partial_r=ap$r,atlas_mean_log2TPM=rowMeans(log2(AT+1)),
  atlas_sd_log2TPM=apply(log2(AT+1),1,sd),atlas_empirical_p=afd$p,
  atlas_pooled_call=(ar>=afd$pos | ar<=-afd$neg))
av[gene_id==FOCAL,atlas_pooled_call:=FALSE]
save(av,'atlas_all_genes.tsv.gz')
logmsg('Atlas thresholds +',afd$pos,'-',afd$neg,'positive fraction',mean(ar[-afi]>0))

# Fixed candidates come solely from the discovery calculation above.
candidate<-tab[stable_discovery==TRUE & gene_id!=FOCAL]
candidate<-merge(candidate,av,by='gene_id',all.x=TRUE)
testids<-candidate[!is.na(atlas_within_r),gene_id];idx<-match(testids,ag)
sgn<-sign(candidate$discovery_within_r[match(testids,candidate$gene_id)])
blocks<-split(seq_len(ncol(AW)),runs$condition);tail<-integer(length(idx))
R<-row_unit(AW[idx,,drop=FALSE]);obs<-ar[idx]*sgn
logmsg('Validating',length(idx),'fixed candidates with',NVALID,'within-condition permutations')
for(start in seq(1,NVALID,by=100)){
 b<-min(100,NVALID-start+1)
 P<-replicate(b,{xx<-AW[afi,];for(j in blocks)xx[j]<-sample(xx[j]);xx/sqrt(sum(xx*xx))})
 rr<-(R%*%P)*sgn;tail<-tail+rowSums(rr>=obs-1e-12)
}
vp<-(tail+1)/(NVALID+1);vq<-p.adjust(vp,'BH')
candidate[,validation_directional_p:=vp[match(gene_id,testids)]]
candidate[,validation_directional_BH:=vq[match(gene_id,testids)]]
candidate[,same_sign:=sign(atlas_within_r)==sign(discovery_within_r)]
candidate[,replicated:=!is.na(validation_directional_BH) & validation_directional_BH<.05 & same_sign]
atlas_loo<-leave_one_min(AW,afi,runs$condition,rep(1,nrow(AW)))
candidate[,atlas_minimum_leave_one_sample_out_r:=atlas_loo[match(gene_id,ag)]]
avloo<-sapply(unique(runs$condition),function(k){j<-which(runs$condition!=k);cr(AW[,j,drop=FALSE],AW[afi,j])})
candidate[,atlas_minimum_leave_one_condition_out_r:=apply(avloo,1,min)[match(gene_id,ag)]]

# Conditional validation is sensitivity, not a new genome-wide discovery screen.
# Freedman-Lane style residual permutations, reprojected on each leave-partner-out score.
conditional_p<-rep(NA_real_,length(idx))
perm_indices<-replicate(NVALID,{ix<-seq_len(ncol(AW));for(bb in blocks)ix[bb]<-sample(bb);ix})
for(k in seq_along(idx)){
 j<-idx[k];s<-ap$score[j,];xr<-ap$x[j,];yr<-ap$y[j,];ss<-sum(s*s)
 P<-matrix(xr[perm_indices],nrow=ncol(AW))
 # yr is already orthogonal to s; only the permuted x norm needs re-projection.
 ps<-as.numeric(crossprod(s,P))
 rr<-as.numeric(crossprod(yr,P))/sqrt(sum(yr*yr)*pmax(sum(xr*xr)-ps*ps/pmax(ss,1e-12),1e-12))
 conditional_p[k]<-(1+sum(rr*sgn[k]>=ap$r[j]*sgn[k]-1e-12))/(NVALID+1)
}
candidate[,conditional_validation_p:=conditional_p[match(gene_id,testids)]]
candidate[,conditional_validation_BH:=p.adjust(conditional_validation_p,'BH')]
candidate[,retains_signature_adjusted_support:=replicated & partial_r_beyond_signature*sign(discovery_within_r)>=.30 &
           atlas_partial_r*sign(discovery_within_r)>0 & conditional_validation_BH<.05]

# Expression/variance-matched background: compare mean signed validation r and percentile ranks.
# Candidates were fixed on discovery data, but this remains retrospective validation.
common<-intersect(tab$gene_id,av$gene_id);bg<-merge(tab[gene_id %in% common],av,by='gene_id')
bg[,expression_bin:=pmin(5L,ceiling(frank(atlas_mean_log2TPM,ties.method='average')/.N*5))]
bg[,variance_bin:=pmin(5L,ceiling(frank(atlas_sd_log2TPM,ties.method='average')/.N*5))]
bg[,bin:=paste(expression_bin,variance_bin,sep='_')]
matched<-list()
for(s in c(1,-1)){
 ids<-candidate[same_sign %in% c(TRUE,FALSE) & sign(discovery_within_r)==s & !is.na(atlas_within_r),gene_id]
 if(!length(ids))next
 tt<-bg[match(ids,gene_id)];pool<-bg[!gene_id %in% c(FOCAL,testids)]
 bins<-split(pool$gene_id,pool$bin)
 draws<-replicate(2000,{sel<-unlist(lapply(split(tt$gene_id,tt$bin),function(ids2){b<-tt$bin[match(ids2[1],tt$gene_id)];sample(bins[[b]],length(ids2),replace=length(bins[[b]])<length(ids2))}));mean(av$atlas_within_r[match(sel,av$gene_id)]*s)})
 ob<-mean(tt$atlas_within_r*s)
 matched[[length(matched)+1]]<-data.table(direction=if(s==1)'positive' else 'negative',n=length(ids),observed_mean_signed_r=ob,
   matched_null_median=median(draws),matched_null_025=quantile(draws,.025),matched_null_975=quantile(draws,.975),
   enrichment_p=(1+sum(draws>=ob))/2001)
}
matched<-rbindlist(matched);save(matched,'validation_matched_background.tsv')
candidate<-merge(candidate,ann[,.(gene_id,v51_symbol,at_symbol,at_gene,at_description,original_robust=robust,original_subcluster=subcluster)],by='gene_id',all.x=TRUE)
candidate[,label:=fifelse(!is.na(at_symbol) & at_symbol!='',at_symbol,sub('Vitvi05_01','',gene_id))]
setorder(candidate,validation_directional_BH,-discovery_within_r)
save(candidate,'fixed_candidates_validation.tsv');save(candidate[replicated==TRUE],'replicated_partners.tsv')
save(candidate[retains_signature_adjusted_support==TRUE],'partners_beyond_signature.tsv')

# Historical overlap check with the originally reported robust-positive set.
old<-merge(ann,av,by='gene_id');n<-old[robust==TRUE & direction=='positive',.N]
hit<-old[robust==TRUE & direction=='positive' & atlas_pooled_call==TRUE & atlas_within_r>0,.N]
K<-old[atlas_pooled_call==TRUE & atlas_within_r>0,.N];N<-nrow(old)
hist<-data.table(method='new_positive_tail_FDR',threshold=afd$pos,original_positive_partners_measured=n,
 atlas_positive_calls=K,overlap=hit,universe=N,expected=n*K/N,fold=hit/(n*K/N),p_hyper=phyper(hit-1,K,N-K,n,lower.tail=FALSE))
# Original summary rounded the absolute-r threshold to 0.84; its grid point 0.845
# recovers the reported 153 genes / 11 overlaps. This is a fixed-cutoff check,
# not a claim that the original random-seed-specific FDR curve was reproduced.
oldK<-old[atlas_within_r>=.845,.N];oldHit<-old[robust==TRUE & direction=='positive' & atlas_within_r>=.845,.N]
hist<-rbind(hist,data.table(method='legacy_fixed_cutoff_check',threshold=.845,original_positive_partners_measured=n,
 atlas_positive_calls=oldK,overlap=oldHit,universe=N,expected=n*oldK/N,fold=oldHit/(n*oldK/N),p_hyper=phyper(oldHit-1,oldK,N-oldK,n,lower.tail=FALSE)))
save(hist,'historical_overlap_check.tsv')

# Pairwise edges among replicated partners, not just focal-to-partner spokes.
# These post-selection edges are descriptive; no edge-level FDR claim is made.
rp<-candidate[replicated==TRUE];rp[,strength:=pmin(abs(discovery_within_r),abs(atlas_within_r))]
setorder(rp,-strength);net<-head(rp,30);ng<-net$gene_id
edges<-data.table()
if(length(ng)>=2){
 dcor<-tcrossprod(row_unit(W[match(ng,genes),,drop=FALSE]));acor<-tcrossprod(row_unit(AW[match(ng,ag),,drop=FALSE]))
 pairs<-which(upper.tri(dcor),arr.ind=TRUE)
 edges<-data.table(gene1=ng[pairs[,1]],gene2=ng[pairs[,2]],r_discovery=dcor[pairs],r_validation=acor[pairs])
 edges[,display_edge:=abs(r_discovery)>=.50 & abs(r_validation)>=.50 & sign(r_discovery)==sign(r_validation)]
}
save(edges,'partner_pairwise_edges.tsv');save(net,'displayed_network_nodes.tsv')
coherence<-list()
for(s in c(1,-1)){
 ids<-rp[sign(discovery_within_r)==s,gene_id];if(length(ids)<2)next
 aa<-match(ids,ag);rr<-tcrossprod(row_unit(AW[aa,,drop=FALSE]));ob<-mean(rr[upper.tri(rr)])
 # Significance not claimed: selected partners share the focal gene by design.
 coherence[[length(coherence)+1]]<-data.table(direction=s,n=length(ids),mean_partner_pairwise_validation_r=ob,
  note='Descriptive after focal-gene selection; shared focal correlation itself induces coherence')
}
save(rbindlist(coherence),'partner_coherence.tsv')

# Record sample-level signature values for plotting and independent auditing.
ss<-data.table(sample_uid=meta$sample_uid[early],group=g,hormone=meta$hormone[early],
 focal_within_z=x,programme_within_z=score)
save(ss,'early_sample_residuals.tsv')
logmsg('Replicated',sum(candidate$replicated),'of',length(testids),'measured fixed candidates; beyond signature',sum(candidate$retains_signature_adjusted_support))
summary<-data.table(metric=c('discovery_genes','stable_discovery','validated_candidates_measured','replicated_candidates','beyond_signature',
 'atlas_genes','atlas_fraction_positive','atlas_median_r','early_signature_r','validation_signature_r'),
 value=c(length(genes)-1,sum(stable),length(testids),sum(candidate$replicated),sum(candidate$retains_signature_adjusted_support),
 length(ag)-1,mean(ar[-afi]>0),median(ar[-afi]),cor(x,score),cor(AW[afi,],colMeans(AW[ag %in% up,,drop=FALSE]))))
save(summary,'summary.tsv');print(summary);print(matched)

# Export publication figures with all partners available in the table, not all labelled.
theme_set(theme_classic(base_size=10)+theme(legend.position='bottom',plot.title=element_text(face='bold')))
plotdata<-merge(tab,av,by='gene_id');plotdata<-plotdata[gene_id!=FOCAL]
plotdata[,status:='Other measured genes']
plotdata[gene_id %in% candidate$gene_id,status:='Stable discovery']
plotdata[gene_id %in% rp$gene_id,status:='Replicated']
cols<-c('Other measured genes'='#bfc4c7','Stable discovery'='#2678ae','Replicated'='#c45236')
p1<-ggplot(plotdata,aes(discovery_within_r,atlas_within_r))+geom_point(aes(colour=status),size=.65,alpha=.5)+
 scale_colour_manual(values=cols)+geom_hline(yintercept=0,colour='grey70',linewidth=.3)+geom_vline(xintercept=0,colour='grey70',linewidth=.3)+
 labs(title='A  Co-expression within replicate groups',x='Discovery: within-group Pearson r (90 cultures)',y='Independent callus: within-condition Pearson r (27 cultures)',colour=NULL)
p2<-ggplot(candidate[!is.na(atlas_partial_r)],aes(partial_r_beyond_signature,atlas_partial_r))+geom_hline(yintercept=0,colour='grey70')+
 geom_vline(xintercept=0,colour='grey70')+geom_point(aes(colour=retains_signature_adjusted_support),size=1.3)+
 scale_colour_manual(values=c('FALSE'='#9ba3a9','TRUE'='#c45236'),labels=c('Other candidates','Supported after signature adjustment'))+
 labs(title='B  After accounting for the positive signature',x='Discovery: partial Pearson r',y='Independent callus: partial Pearson r',colour=NULL)
lab<-head(rp,8)
p1<-p1+geom_text_repel(data=lab,aes(discovery_within_r,atlas_within_r,label=label),inherit.aes=FALSE,size=2.4,max.overlaps=20)
p3<-ggplot(ss,aes(focal_within_z,programme_within_z,colour=hormone))+geom_point(size=1.8,alpha=.8)+
 labs(title=sprintf('C  Early treated cultures: signature association (r = %.2f)',cor(x,score)),
 x='VvECPG1: deviation from group mean (z scale)',y='Positive signature: deviation from group mean',colour=NULL)
p4<-ggplot(candidate[!is.na(atlas_within_r)],aes(atlas_within_r,fill=factor(sign(discovery_within_r))))+geom_histogram(bins=25,alpha=.6,position='identity')+
 geom_vline(xintercept=median(ar[-afi]),linetype=2)+scale_fill_manual(values=c('-1'='#2678ae','1'='#c45236'),labels=c('Negative discovery partners','Positive discovery partners'))+
 labs(title='D  Independent correlations of fixed candidates',x='Independent within-condition Pearson r',y='Number of discovery candidates',fill=NULL)
fig<-(p1+p2)/(p3+p4)+plot_annotation(caption='Exploratory follow-up. Candidate discovery uses 6–24 h treated cultures only; validation uses PRJNA732451.\nDashed line in D: transcriptome-wide median. Signature adjustment excludes the focal gene and each tested partner.')
for(ext in c('pdf','png'))ggsave(file.path(OUT,paste0('coexpression_validation.',ext)),fig,width=13,height=10,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else 'png')

# A compact, descriptive network from edges reproduced in both datasets.
if(nrow(net)>0){
 ids<-c(FOCAL,ng);labels<-c('VvECPG1',net$label);n<-length(ng)
 coords<-data.table(gene_id=ids,label=labels,x=c(0,cos(seq(0,2*pi,length.out=n+1)[seq_len(n)])),
 y=c(0,sin(seq(0,2*pi,length.out=n+1)[seq_len(n)])))
 spoke<-data.table(gene1=FOCAL,gene2=ng,type='Focal association')
 others<-data.table()
 if(nrow(edges)){
  ed<-copy(edges[display_edge==TRUE]);ed[,strength:=pmin(abs(r_discovery),abs(r_validation))]
  # Limit clutter without changing the complete exported pairwise-edge table.
  directed<-rbind(ed[,.(gene1,gene2,strength)],ed[,.(gene1=gene2,gene2=gene1,strength)])
  setorder(directed,gene1,-strength);dd<-directed[,head(.SD,2),by=gene1]
  selected<-unique(dd[,.(a=pmin(gene1,gene2),b=pmax(gene1,gene2))])
  others<-selected[,.(gene1=a,gene2=b,type='Partner association')]
 }
 ee<-rbind(spoke,others,fill=TRUE);ee<-merge(ee,coords[,.(gene1=gene_id,x,y)],by='gene1');ee<-merge(ee,coords[,.(gene2=gene_id,xend=x,yend=y)],by='gene2')
 pn<-ggplot()+geom_segment(data=ee,aes(x,y,xend=xend,yend=yend,alpha=type),colour='#9ca4aa')+
  scale_alpha_manual(values=c('Focal association'=.5,'Partner association'=.18))+geom_point(data=coords,aes(x,y),size=3,colour='#c45236')+
  geom_label_repel(data=coords,aes(x,y,label=label),size=2.4,max.overlaps=Inf,seed=1)+coord_equal()+theme_void()+theme(legend.position='none')+
  labs(title='VvECPG1 and reproducible co-expression partners',subtitle='Up to 30 partners with strongest minimum correlation across datasets',
  caption='Partner edges: up to two strongest links per node, |r| >= 0.50 with the same sign in both datasets.\nLabels use Arabidopsis homolog symbols where available. Edges are descriptive, not direct regulation.')
 for(ext in c('pdf','png'))ggsave(file.path(OUT,paste0('replicated_neighbourhood.',ext)),pn,width=10,height=10,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else 'png')
}
saveRDS(list(discovery_genes=genes,discovery_within=W,discovery_groups=g,validation_genes=ag,
 validation_within=AW,validation_groups=runs$condition,parameters=list(seed=20261108,NPERM=NPERM,NVALID=NVALID,NBOOT=NBOOT)),file.path(OUT,'analysis_matrices.rds'))
capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
logmsg('DONE. Outputs:',OUT)
