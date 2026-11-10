#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
ROOT <- normalizePath('.')
SRC <- file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "functional")
OUT <- file.path(ROOT,'results/figures/functional')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
d <- fread(file.path(SRC,'displayed_GO_enrichment.tsv'))
jj <- fread(file.path(SRC,'displayed_ranked_processes.tsv'))
laborder <- c('DNA damage response','DNA repair','DNA metabolic process')
d[,term_name:=factor(term_name,levels=rev(laborder))]
d[,background:=factor(background,levels=c('Shared measurable genes','Positive discovery candidates'))]
theme_set(theme_classic(base_size=10)+theme(legend.position='bottom',
 plot.title=element_text(face='bold',size=11,margin=margin(b=4)),
 plot.subtitle=element_text(size=9,margin=margin(b=4)),
 axis.title=element_text(size=10),legend.text=element_text(size=9),
 legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),
 legend.key.height=unit(10,'pt'),legend.spacing.y=unit(0,'pt'),
 plot.margin=margin(3,5,3,3),panel.spacing=unit(7,'pt')))
p1 <- ggplot(d,aes(fold_enrichment,term_name,colour=significance))+
 geom_point(size=2.7)+geom_text(aes(label=paste0('  n=',k)),hjust=0,size=2.8,show.legend=FALSE)+
 scale_colour_manual(values=c('BH < 0.05'='#b94a31','BH >= 0.05'='#7d8992'))+
 scale_x_continuous(limits=c(0,10),breaks=seq(0,10,2))+facet_wrap(~background)+
 labs(title='A  DNA-maintenance annotations among 106 replicated partners',
 x='Fold enrichment among GO-annotated genes',y=NULL,colour=NULL)
rd <- rbind(jj[,.(term_name,term_id,dataset='Discovery',AUC=AUC_discovery)],
 jj[,.(term_name,term_id,dataset='Independent callus',AUC=AUC_validation)])
rd[,term_name:=factor(term_name,levels=rev(laborder))]
p2 <- ggplot(rd,aes(AUC,term_name,colour=dataset))+
 geom_vline(xintercept=.5,linetype=2,colour='grey65')+
 geom_point(size=2.7,position=position_dodge(width=.35))+
 scale_colour_manual(values=c('Discovery'='#2678ae','Independent callus'='#c45236'))+
 scale_x_continuous(limits=c(.45,.85),breaks=seq(.5,.8,.1))+
 labs(title='B  Rank-based effects for the same categories',
 subtitle='No GO BP or KEGG set passed the corrected cross-dataset conjunction test',
 x='AUC (0.5 = no rank advantage)',y=NULL,colour=NULL)
fig <- p1/p2+plot_layout(heights=c(1,1))+plot_annotation(
 caption='GO terms overlap; 55 of 106 partners have transferred GO BP annotations.\nArabidopsis best-hit annotations do not establish grapevine protein function. Categories in B were selected from A.',
 theme=theme(plot.caption=element_text(size=8,hjust=0,margin=margin(t=3)),plot.margin=margin(4,4,4,4)))
for(ext in c('pdf','png'))ggsave(file.path(OUT,paste0('functional_panels.',ext)),fig,width=8.5,height=5,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else 'png')
fwrite(d,file.path(OUT,'verified_enrichment_points.tsv'),sep='\t')
fwrite(rd,file.path(OUT,'verified_ranked_points.tsv'),sep='\t')
cat('Redrew',nrow(d),'enrichment points and',nrow(rd),'ranked points; original estimates and thresholds unchanged.\n')
