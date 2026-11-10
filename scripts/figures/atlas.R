suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
R <- normalizePath('.');Q <- file.path('results','figures','atlas');dir.create(Q,recursive=TRUE,showWarnings=FALSE)
s <- fread(file.path(R,'inputs/atlas/focal_atlas_per_sample.tsv')); set.seed(1)
order <- c('embryogenic callus','unlabelled calli (inferred embryogenic)','undifferentiated callus, 40 d (Sangiovese)','non-embryogenic callus','callus / cell culture','somatic embryo','seed / zygotic embryo','flower / inflorescence','berry / fruit','berry skin','berry pulp / flesh','bud','shoot apex','shoot / seedling / whole plant','leaf','tendril','stem / cane / xylem / phloem','root')
group <- setNames(c(rep('Culture',6),rep('Reproductive',5),rep('Vegetative',7)),order)
co <- c('Culture'='#2a78d6','Reproductive'='#eb6834','Vegetative'='#1baf7a');refs <- c('UBQ2'='#6250d6','PP2A-A'='#eda100','GAPDH'='#e87ba4')
s[,group:=factor(group[tissue_class],levels=names(co))]
short <- setNames(order,order);short['callus / cell culture'] <- 'other callus / cell culture';short['undifferentiated callus, 40 d (Sangiovee)'] <- 'undifferentiated callus, 40 d †'
short['undifferentiated callus, 40 d (Sangiovese)'] <- 'undifferentiated callus, 40 d †';short['unlabelled calli (inferred embryogenic)'] <- 'unlabelled calli (inferred EC) *'
export <- function(fig,name,w,h){
 for(ext in c('png','pdf','svg')){
  f <- file.path(Q,paste0(name,'.',ext))
  ggsave(f,fig,width=w,height=h,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else if(ext=='svg')svg else 'png')
 }
}
primary <- s[tissue_class!='unlabelled calli (inferred embryogenic)'];ord <- order[order %in% primary$tissue_class]
primary[,ix:=match(tissue_class,ord)]
b <- primary[,.(v=median(focal_tpm)),by=.(tissue_class,study_accession,ix)]
t <- b[,.(v=median(v)),by=.(tissue_class,ix)]
n <- primary[,.(samples=.N,studies=uniqueN(study_accession)),by=.(tissue_class,ix)]
band <- data.table(ix=seq_along(ord))[ix%%2==0]
groups <- primary[,.(lo=min(ix)-.4,hi=max(ix)+.4,mid=mean(range(ix))),by=group]
main <- ggplot()+geom_rect(data=band,aes(xmin=ix-.5,xmax=ix+.5,ymin=.38,ymax=110000),fill='#ecebe6')+
 geom_point(data=primary,aes(ix,focal_tpm+.5,colour=group),position=position_jitter(width=.23,height=0,seed=1),size=1.05,alpha=.5)+
 geom_segment(data=b,aes(x=ix-.21,xend=ix+.21,y=v+.5,yend=v+.5),linewidth=.5)+
 geom_point(data=b,aes(ix,v+.5,shape='study median'),size=3,colour='black')+
 geom_point(data=t,aes(ix,v+.5,shape='median of study medians'),fill='white',size=2.5,stroke=.9)+
 geom_vline(xintercept=c(5.5,10.5),colour='#52514e',linewidth=.3)+
 geom_text(data=n,aes(ix,160000,label=paste0(samples,'\n',studies)),size=3.2,colour='#52514e',lineheight=.9)+
 geom_segment(data=groups,aes(x=lo,xend=hi,y=450000,yend=450000,colour=group),linewidth=.8,show.legend=FALSE)+
 geom_text(data=groups,aes(mid,700000,label=group,colour=group),fontface='bold',size=4,show.legend=FALSE)+
 annotate('text',x=.3,y=260000,label='samples /\nstudies',hjust=1,size=3.1,fontface='bold',lineheight=.9)+
 scale_colour_manual(values=co,name=NULL)+scale_shape_manual(values=c('study median'=95,'median of study medians'=23),name=NULL)+scale_y_log10(breaks=10^(0:5),labels=c('1','10','100','1000','10000','100000'),expand=c(0,0))+
 scale_x_continuous(breaks=seq_along(ord),labels=unname(short[ord]),expand=expansion(add=.6))+
 coord_cartesian(ylim=c(.38,110000),clip='off')+
 labs(title='VvECPG1 expression across public grapevine RNA-seq data',x=NULL,y='TPM (log scale; 0 at bottom edge)')+
 theme_classic(base_size=11)+theme(panel.background=element_rect(fill='#f7f7f4',colour=NA),axis.line=element_blank(),
 panel.grid.major.y=element_line(colour='#dcdbd5',linewidth=.3),axis.ticks.x=element_blank(),
 axis.text.x=element_text(angle=55,hjust=1),legend.position='top',legend.margin=margin(0,0,0,0),
 plot.title=element_text(face='bold',size=15,margin=margin(b=10)),plot.margin=margin(15,8,35,12),
 legend.box.margin=margin(b=75))
# Put header counts within a reserved upper margin, separate from the log-scaled data.
export(main,'focal_atlas_main',10,8.5)
# Horizontal atlas panels retain BioSample dots, study medians and study-weighted summaries.
horizontal <- function(dat,full=FALSE){
 ord <- order[order %in% dat$tissue_class]; dat <- copy(dat);dat[,iy:=length(ord)+1-match(tissue_class,ord)]
 bands <- data.table(iy=seq_along(ord))[iy%%2==0]
 base <- function(log=FALSE,low=1) ggplot()+geom_rect(data=bands,aes(ymin=iy-.5,ymax=iy+.5,xmin=stage(start=low,after_scale=-Inf),xmax=stage(start=low,after_scale=Inf)),fill='#ecebe6')+
  scale_y_continuous(breaks=seq_along(ord),labels=unname(short[rev(ord)]),expand=expansion(add=.6))+
  theme_classic(base_size=10)+theme(axis.line=element_blank(),axis.ticks.y=element_blank(),panel.background=element_rect(fill='#f7f7f4',colour=NA),panel.grid.major.x=element_line(colour='#dcdbd5',linewidth=.3),plot.title=element_text(face='bold',size=11),plot.margin=margin(6,8,6,6))+
  labs(y=NULL)
 strip <- function(col,log=FALSE){
  z <- copy(dat);z[,v:=get(col)+if(log).5 else 0];m <- z[,.(v=median(v)),by=.(iy,study_accession)];tt<-m[,.(v=median(v)),by=iy]
  pp<-base(log,min(z$v))+geom_point(data=z,aes(v,iy,colour=group),position=position_jitter(width=0,height=.18,seed=1),size=1.2,alpha=.5)+
   geom_segment(data=m,aes(x=v,xend=v,y=iy-.2,yend=iy+.2),linewidth=.5)+geom_point(data=tt,aes(v,iy),shape=23,fill='white',size=2.3,stroke=.8)+scale_colour_manual(values=co,name=NULL)+theme(legend.position='bottom')
  if(log)pp<-pp+scale_x_log10() else pp<-pp+geom_vline(xintercept=0,linetype=2,colour='grey50')
  pp
 }
 a<-strip('focal_tpm',TRUE)+labs(title='A. Focal gene',x='TPM (log scale; 0 at left edge)')
 b<-strip('log2_focal_over_UBQ2')+labs(title=paste0(if(full)'B' else 'A','. Focal gene relative to UBQ2'),subtitle='VIT_19s0177g00070',x='log2(focal / UBQ2)')
 rr<-melt(dat,id.vars=c('iy','tissue_class','study_accession','sample_accession'),measure.vars=paste0(names(refs),'_tpm'),variable.name='reference',value.name='v');rr[,reference:=sub('_tpm$','',reference)];rr[,off:=c('UBQ2'=.22,'PP2A-A'=0,'GAPDH'=-.22)[reference]]
 rm<-rr[,.(v=median(v)),by=.(iy,study_accession,reference,off)][,.(v=median(v)),by=.(iy,reference,off)]
 c<-base(TRUE,min(rr$v+.5))+geom_point(data=rr,aes(v+.5,iy+off,colour=reference),position=position_jitter(width=0,height=.05,seed=1),alpha=.35,size=.9)+
  geom_point(data=rm,aes(v+.5,iy+off,colour=reference),size=2.5,stroke=.5)+scale_colour_manual(values=refs,name=NULL)+scale_x_log10()+
  theme(axis.text.y=element_blank(),legend.position='bottom')+labs(title=paste0(if(full)'C' else 'B','. Reference genes'),x='TPM (log scale)')
 n<-dat[,.(samples=.N,studies=uniqueN(study_accession)),by=iy]
 c<-c+geom_text(data=n,aes(x=Inf,y=iy,label=paste0(samples,' / ',studies)),hjust=-.12,inherit.aes=FALSE,size=2.5)+coord_cartesian(clip='off')+theme(plot.margin=margin(6,65,6,6))
 if(full){b<-b+theme(axis.text.y=element_blank()); fig<-a+b+c+plot_layout(widths=c(1.2,1,1))}else fig<-b+c+plot_layout(widths=c(1,1))
 fig+plot_annotation(caption='Dots: BioSamples (technical runs merged). Black ticks: study medians; diamonds: median of study medians.\nBetween-study comparisons are descriptive. * Inferred embryogenic identity; † Sangiovese 40-d undifferentiated callus.',theme=theme(plot.caption=element_text(size=8,hjust=0)))
}
export(horizontal(s,TRUE),'focal_atlas_full',17,9)
export(horizontal(primary,FALSE),'focal_atlas_supp',12.5,9)
cat('Rebuilt primary atlas:',nrow(primary),'BioSamples;',sum(primary$n_runs),'runs.\n')
