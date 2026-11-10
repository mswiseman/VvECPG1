library(data.table)
library(ggplot2)
library(ggrepel)
ROOT <- normalizePath('.')
OUT <- file.path(ROOT,'results/figures/network')
SRC <- file.path(Sys.getenv("VVECPG1_RESULTS", "results"), "coexpression")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
net <- fread(file.path(SRC,'displayed_network_nodes.tsv'))
edges <- fread(file.path(SRC,'partner_pairwise_edges.tsv'))
FOCAL <- 'Vitvi05_01chr01g05750'; ng <- net$gene_id; n <- length(ng)
angle <- seq(0,2*pi,length.out=n+1)[seq_len(n)]
# Fixed plotting bounds make the reduced radius visibly shorten edges.
radius <- 0.66
coords <- data.table(gene_id=c(FOCAL,ng),label=c('VvECPG1',net$label),
 x=c(0,radius*cos(angle)),y=c(0,radius*sin(angle)))
spoke <- data.table(gene1=FOCAL,gene2=ng,type='Focal association')
ed <- copy(edges[display_edge==TRUE]);ed[,strength:=pmin(abs(r_discovery),abs(r_validation))]
directed <- rbind(ed[,.(gene1,gene2,strength)],ed[,.(gene1=gene2,gene2=gene1,strength)])
setorder(directed,gene1,-strength); dd <- directed[,head(.SD,2),by=gene1]
selected <- unique(dd[,.(a=pmin(gene1,gene2),b=pmax(gene1,gene2))])
others <- selected[,.(gene1=a,gene2=b,type='Partner association')]
ee <- rbind(spoke,others,fill=TRUE)
ee <- merge(ee,coords[,.(gene1=gene_id,x,y)],by='gene1')
ee <- merge(ee,coords[,.(gene2=gene_id,xend=x,yend=y)],by='gene2')
pn <- ggplot()+geom_segment(data=ee,aes(x,y,xend=xend,yend=yend,alpha=type),colour='#9ca4aa')+
 scale_alpha_manual(values=c('Focal association'=.5,'Partner association'=.18))+
 geom_point(data=coords,aes(x,y),size=3,colour='#c45236')+
 geom_label_repel(data=coords,aes(x,y,label=label),size=2.8,max.overlaps=Inf,seed=1,
  box.padding=.28,point.padding=.15,min.segment.length=0)+
 coord_equal(xlim=c(-1.02,1.02),ylim=c(-.94,.94),expand=FALSE,clip='off')+
 theme_void()+theme(legend.position='none',plot.title=element_text(size=14),
 plot.subtitle=element_text(size=10),plot.caption=element_text(size=8),plot.margin=margin(12,16,12,16))+
 labs(title='VvECPG1 and reproducible co-expression partners',
 subtitle='Up to 30 partners with strongest minimum correlation across datasets',
 caption='Partner edges: up to two strongest links per node, |r| >= 0.50 with the same sign in both datasets.\nLabels use Arabidopsis homolog symbols where available. Edges are descriptive, not direct regulation.')
for(ext in c('pdf','png'))ggsave(file.path(OUT,paste0('network_panel.',ext)),pn,width=8,height=8,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else 'png')
fwrite(ee,file.path(OUT,'displayed_edges_verification.tsv'),sep='\t')
cat('Preserved',n,'partner nodes;',nrow(spoke),'focal edges;',nrow(others),'partner edges.\n')
