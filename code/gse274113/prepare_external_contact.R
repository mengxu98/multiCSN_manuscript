suppressPackageStartupMessages({library(data.table);library(rtracklayer)})
root <- "data/gse274113/raw_external_downloads"
out <- "data/gse274113/external_resources"
chain_file <- file.path(root,"hg19ToHg38.over.chain")
if(!file.exists(chain_file)) {
 con <- gzfile(paste0(chain_file,".gz"),"rt")
 writeLines(readLines(con),chain_file);close(con)
}
chain <- import.chain(chain_file)
d <- fread(file.path(root,"HSPC_loops_hg19.tsv"))
mapped <- lapply(1:2,function(k) {
 prefix<-paste0("anchor",k,"_")
 gr<-GRanges(d[[paste0(prefix,"chrom")]],IRanges(d[[paste0(prefix,"start")]]+1L,d[[paste0(prefix,"end")]]))
 liftOver(gr,chain)
})
keep <- lengths(mapped[[1]])==1L & lengths(mapped[[2]])==1L
result <- d[keep]
for(k in 1:2){
 gr<-unlist(mapped[[k]][keep],use.names=FALSE)
 set(result,j=paste0("anchor",k,"_chrom"),value=as.character(seqnames(gr)))
 set(result,j=paste0("anchor",k,"_start"),value=start(gr))
 set(result,j=paste0("anchor",k,"_end"),value=end(gr))
}
result[,coordinate_system:="1-based closed"]
result[,build:="hg38"]
fwrite(result,file.path(out,"HSPC_loops_hg38.tsv.gz"),sep="\t")
d[,`:=`(anchor1_mappings=lengths(mapped[[1]]),anchor2_mappings=lengths(mapped[[2]]),retained=keep)]
fwrite(d,file.path(root,"contact_liftover_accounting.tsv"),sep="\t")
print(d[,.(input_loops=.N,retained=sum(retained)),by=source_sheet])
