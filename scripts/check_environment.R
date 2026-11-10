#!/usr/bin/env Rscript
packages <- c('data.table','Matrix','ggplot2','patchwork','ggrepel')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
cat(R.version.string,'\n')
for(p in packages)cat(p,as.character(packageVersion(p)),'\n')
if(!nzchar(Sys.which('gzip')))stop('gzip executable is required')
if(!capabilities('cairo'))stop('Cairo graphics support is required for PDF export')
