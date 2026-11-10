focal<-'Vitvi05_01chr01g05750';partner<-'Vitvi05_01chr01g05760';out<-'results/paralog';dir.create(out,recursive=TRUE,showWarnings=FALSE)
a<-readRDS('inputs/focal_network/context_C_adjusted.rds');m<-a$adj
x<-as.numeric(m[focal,]);y<-as.numeric(m[partner,])
bw<-function(z){u<-(z-median(z))/(9*mad(z,constant=1));w<-ifelse(abs(u)<1,(1-u*u)^2,0);v<-(z-median(z))*w;v/sqrt(sum(v*v))}
r<-sum(bw(x)*bw(y));cat('Legacy bicor direct standard formula:',format(r,digits=16),'\n')
cat('Max outlier fractions',mean(abs((x-median(x))/(9*mad(x,constant=1)))>=1),mean(abs((y-median(y))/(9*mad(y,constant=1)))>=1),'\n')
write.table(data.frame(sample_uid=colnames(m),focal_adjusted=x,partner_adjusted=y),file.path(out,'legacy_adjusted_pair.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
b<-readRDS(file.path(file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression"),"analysis_matrices.rds"))
w<-b$discovery_within;v<-b$validation_within
r1<-cor(w[focal,],w[partner,]);r2<-cor(v[focal,],v[partner,]);cat('New discovery Pearson',format(r1,digits=16),'n',ncol(w),'\n');cat('Independent validation Pearson',format(r2,digits=16),'n',ncol(v),'\n')
write.table(data.frame(sample_uid=colnames(w),group=b$discovery_groups,focal=w[focal,],partner=w[partner,]),file.path(out,'current_discovery_pair.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(sample_uid=colnames(v),group=b$validation_groups,focal=v[focal,],partner=v[partner,]),file.path(out,'current_validation_pair.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(context=c('legacy_126_libraries','current_90_treated_discovery','independent_27_samples'),method=c('biweight_midcorrelation_standard_formula','within_group_Pearson','within_condition_Pearson'),n=c(length(x),ncol(w),ncol(v)),correlation=c(r,r1,r2)),file.path(out,'correlation_verification.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
