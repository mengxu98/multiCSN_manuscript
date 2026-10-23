suppressPackageStartupMessages(library(data.table))
source("code/functions/go_term_selection.R")
source("code/requirements/prepare_go_resources.R")

export_go_annotation <- function(destination = "results/gse194122/go_annotation", resource_root = "data/resources/fig5_go_3.22.0") {
prepare_display_go_similarity(resource_root)
if (dir.exists(destination)) stop("Output already exists: ", destination)
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
resource_path <- file.path(resource_root, "cache/af3b10e4c2379820a80657c5a6685dc4.Rcache")
resource <- R.cache::loadCache(pathname = resource_path)
fwrite(unique(as.data.table(resource$TERM2GENE)[, .(go_id = Term, target = symbol)]), file.path(destination, "TERM2GENE.tsv.gz"), sep = "\t")
fwrite(unique(as.data.table(resource$TERM2NAME)[, .(go_id = Term, term = Name)]), file.path(destination, "TERM2NAME.tsv"), sep = "\t")
similarity_path <- file.path(resource_root, "display_wang_similarity.rds")
go_similarity <- readRDS(similarity_path)
x <- fread("results/gse194122/donor_stability/stable_rewiring_GO_BP.tsv.gz")
a <- semantic_select(copy(x[regulator == "HMGA1" & transition_order == 4 & topology_class == "lost" & p_adjust < 0.05]), go_similarity)
b <- merge(x[regulator == "BNC2" & transition_order == 3 & topology_class == "gained" & p_adjust < 0.05], x[regulator == "BNC2" & transition_order == 4 & topology_class == "lost" & p_adjust < 0.05], by = "go_id", suffixes = c("_gain", "_loss"))
b[, worst_q := pmax(p_adjust_gain, p_adjust_loss)]
b <- semantic_select(b, go_similarity, order_col = "worst_q", overlap_col = "overlap_genes_gain")
y <- rbindlist(list(a[,.(regulator,transition_order,topology_class,go_id,term)],b[,.(regulator="BNC2",transition_order=3L,topology_class="gained",go_id,term=term_gain)],b[,.(regulator="BNC2",transition_order=4L,topology_class="lost",go_id,term=term_loss)]))
fwrite(y,file.path(destination, "displayed_primary_terms.tsv"),sep="\t")
print(y)

metadata <- list(version = as.character(resource$version), ontology = "BP",
  species = "Homo_sapiens", similarity = "Wang", primary_display_terms = "displayed_primary_terms.tsv")
jsonlite::write_json(metadata, file.path(destination, "metadata.json"), auto_unbox = TRUE, pretty = TRUE)
}

if (sys.nframe() == 0L) export_go_annotation()
