reevaluate_networks <- function(result_root = "results/beeline") {
# Re-evaluate all 11 methods from networks using a common metric implementation.
# Evaluate networks on any complete execution tree.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
source('code/functions/functions.R')
base <- normalizePath(getwd(), mustWork=TRUE)
result_root <- normalizePath(result_root, mustWork = TRUE)
out <- file.path(result_root, 'evaluation')
if (dir.exists(out)) stop('Output already exists: ', out)
truths <- unlist(lapply(c('Curated','Synthetic'), function(x) list.files(file.path(base,'data/beeline',x), pattern='^refNetwork.csv$', recursive=TRUE, full.names=TRUE)))
methods <- c('multiCSN','GENIE3','GNIPLR','GRNBoost2','LEAP','PIDC','PPCOR','SCODE','SCORPION','SCSGL','SINCERITIES')
rows <- list(); i <- 0L
for (truth_path in truths) {
 case <- basename(dirname(truth_path)); group <- basename(dirname(dirname(truth_path)))
 typ <- tolower(basename(dirname(dirname(dirname(truth_path)))))
 gt <- as.data.frame(fread(truth_path))
 for (method in methods) {
  path <- file.path(result_root,'networks',paste0(method,'__',case,'.csv'))
   nt <- fread(path)
  metric <- fast_standard_metric_summary(nt,gt)
  i <- i+1L
  rows[[i]] <- c(list(Case=case,DatasetType=typ,DatasetGroup=group,Method=method),metric)
 }
}
dir.create(out, recursive=TRUE)
fwrite(rbindlist(rows),file.path(out,'case_metrics.tsv'),sep='\t')
}



if (sys.nframe() == 0L) reevaluate_networks()
