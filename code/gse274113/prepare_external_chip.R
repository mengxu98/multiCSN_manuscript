suppressPackageStartupMessages({library(data.table);library(rtracklayer)})
root <- "data/gse274113/raw_external_downloads"
out <- "data/gse274113/external_resources"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
chain_file <- file.path(root,"hg18ToHg38.over.chain")
if(!file.exists(chain_file)) {
 con <- gzfile(paste0(chain_file,".gz"),"rt")
 writeLines(readLines(con),chain_file);close(con)
}
chain <- import.chain(chain_file)
summaries <- list(); tables <- list()
for(tf in c("GATA1","NFE2")) {
 raw <- fread(file.path(root,paste0(tf,"_hg18.bed.gz")),header=FALSE)
 gr <- GRanges(raw[[1]],IRanges(raw[[2]]+1L,raw[[3]]))
 mapped <- liftOver(gr,chain)
 keep <- lengths(mapped)==1L
 gr38 <- unlist(mapped[keep],use.names=FALSE)
 tables[[tf]] <- data.table(regulator=tf,chrom=as.character(seqnames(gr38)),start=start(gr38),end=end(gr38),source_row=which(keep),source_build="hg18",build="hg38",coordinate_system="1-based closed")
 summaries[[tf]] <- data.table(regulator=tf,input_rows=nrow(raw),unique_mappings=sum(keep),unmapped=sum(lengths(mapped)==0),multiple_mappings=sum(lengths(mapped)>1))
}
fwrite(rbindlist(tables),file.path(out,"primary_erythroid_ChIP_peaks_hg38.tsv.gz"),sep="\t")
fwrite(rbindlist(summaries),file.path(root,"chip_liftover_summary.tsv"),sep="\t")
print(rbindlist(summaries))
