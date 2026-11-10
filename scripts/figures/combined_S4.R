#!/usr/bin/env Rscript
# Reassemble the three panels; layout may differ from the manually composed slide.
source('scripts/figures/network.R')
source('scripts/figures/functional.R')
dir.create('results/figures/combined',recursive=TRUE,showWarnings=FALSE)
pn <- pn + labs(title='A. VvECPG1 and reproducible co-expression partners',subtitle=NULL,caption=NULL)
p1 <- p1 + labs(title='B. DNA-maintenance annotations among 106 replicated partners',
               x='Fold enrichment among GO-annotated replicated partners')
p2 <- p2 + labs(title='C. Rank-based effects for the same categories',subtitle=NULL)
combined <- (pn | (p1 / p2)) + patchwork::plot_layout(widths=c(1,1.25)) +
 patchwork::plot_annotation(caption=paste(
 'GO terms overlap; 55 of 106 partners have transferred GO BP annotations. Categories in C were selected from B.',
 'No GO BP or KEGG set passed the corrected cross-dataset conjunction test.',
 'Homology-transferred annotations and descriptive network edges do not establish grapevine protein function or direct regulation.',sep='\n'),
 theme=theme(plot.caption=element_text(size=8,hjust=0)))
for(ext in c('pdf','png')) ggsave(paste0('results/figures/combined/Figure_S4.',ext),combined,
 width=16,height=8,dpi=300,bg='white',device=if(ext=='pdf')cairo_pdf else 'png')
