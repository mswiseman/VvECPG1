#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
O<-file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "functional")
readtab<-function(f)if(grepl('gz$',f))fread(cmd=paste('gzip -dc',shQuote(f)))else fread(f)
a<-readtab(file.path(O,'enrichment_with_homolog_sensitivity.tsv.gz'))
j<-readtab(file.path(O,'ranked_gene_sets_conjunction.tsv'))
rep<-readtab(file.path(file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression"),"replicated_partners.tsv"))
ann<-readtab(file.path(O,'transferred_gene_annotations.tsv.gz'))
tair_names<-readtab('inputs/functional_annotations/tair_current_gene_names.tsv')
ids<-a[comparison=='replicated_vs_common' & ontology=='BP' & unit=='gene' & BH<.05,term_id]
laborder<-c('DNA damage response','DNA repair','DNA metabolic process')
d<-a[comparison %in% c('replicated_vs_common','replicated_vs_positive_candidates') & ontology=='BP' & unit=='gene' & term_id %in% ids]
d[,term_name:=factor(term_name,levels=rev(laborder))]
d[,background:=factor(ifelse(comparison=='replicated_vs_common','Shared measurable genes','Positive discovery candidates'),levels=c('Shared measurable genes','Positive discovery candidates'))]
d[,significance:=ifelse(BH<.05,'BH < 0.05','BH >= 0.05')]
theme_set(theme_classic(base_size=11)+theme(legend.position='bottom',plot.title=element_text(face='bold')))
p1<-ggplot(d,aes(fold_enrichment,term_name,colour=significance))+geom_point(size=3)+
 geom_text(aes(label=paste0('  n=',k)),hjust=0,size=3,show.legend=FALSE)+
 scale_colour_manual(values=c('BH < 0.05'='#b94a31','BH >= 0.05'='#7d8992'))+
 scale_x_continuous(limits=c(0,10),breaks=seq(0,10,2))+facet_wrap(~background)+
 labs(title='A  DNA-maintenance annotations among 106 replicated partners',x='Fold enrichment among GO-annotated genes',y=NULL,colour=NULL)
jj<-j[ontology=='BP' & term_id %in% ids]
rd<-rbind(jj[,.(term_name,term_id,dataset='Discovery',AUC=AUC_discovery)],
          jj[,.(term_name,term_id,dataset='Independent callus',AUC=AUC_validation)])
rd[,term_name:=factor(term_name,levels=rev(laborder))]
p2<-ggplot(rd,aes(AUC,term_name,colour=dataset))+geom_vline(xintercept=.5,linetype=2,colour='grey65')+
 geom_point(size=3,position=position_dodge(width=.35))+scale_colour_manual(values=c('Discovery'='#2678ae','Independent callus'='#c45236'))+
 scale_x_continuous(limits=c(.45,.85),breaks=seq(.5,.8,.1))+
 labs(title='B  Full-ranking effects for the same categories',subtitle='No GO or KEGG process passed the rank-based conjunction test at BH < 0.05',
 x='Probability a process gene ranks above a background gene',y=NULL,colour=NULL)
fig<-p1/p2+plot_annotation(caption='GO terms overlap and are not independent findings. Only 55 of 106 partners carry a transferred GO BP annotation.\nAnnotations use Arabidopsis best-hit homology; arrows of regulation and VvECPG1 function are not inferred. Panel B categories were selected from panel A.')
for(ext in c('pdf','png'))ggsave(file.path(O,paste0('functional_summary.',ext)),fig,width=12,height=7.5,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else 'png')
fwrite(d,file.path(O,'displayed_GO_enrichment.tsv'),sep='\t')
fwrite(jj,file.path(O,'displayed_ranked_processes.tsv'),sep='\t')
support<-ann[gene_id %in% rep$gene_id & term_id %in% ids]
support<-unique(support[,.(gene_id,term_id)])
support[,present:=TRUE]
wide<-dcast(support,gene_id~term_id,value.var='present',fill=FALSE)
gt<-merge(rep[,.(gene_id,at_gene,original_homolog_symbol=label,in_positive_signature,discovery_within_r,
 atlas_within_r,validation_directional_BH)],wide,by='gene_id')
gt<-merge(gt,tair_names,by='at_gene',all.x=TRUE)
setorder(gt,tair_symbol)
fwrite(gt,file.path(O,'DNA_maintenance_supporting_genes.tsv'),sep='\t')
fwrite(a[unit=='gene',.(comparison,ontology,term_id,term_name,k,n_annotated,K,N_annotated,fold_enrichment,BH,homolog_BH,
 supporting_grape_genes,supporting_homologs)],file.path(O,'Supplementary_Table_functional_enrichment.tsv.gz'),sep='\t')
