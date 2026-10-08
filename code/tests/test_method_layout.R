method <- "multiCSN"
source("code/beeline/unified_worker.R")
stopifnot(identical(method, "multiCSN"), setequal(names(method_runners),
  c("GENIE3", "SCORPION", "PIDC", "LEAP", "PPCOR", "SCODE", "SINCERITIES")),
  is.function(run_multiCSN))

check_vendor_loading <- function() {
  previous <- getwd()
  on.exit(setwd(previous), add = TRUE)
  setwd(tempdir())
  stopifnot(file.exists(find_scode_script()),
    is.function(sincerities_cmtest2()), is.function(sincerities_ranked_predictions()),
    ensure_pidc_available() %in% c("vendor", "package"))
}
check_vendor_loading()
# Exercise the two adapters that convert a fitted matrix into edge tables.
# Loading alone would not catch a missing shared converter dependency.
expression <- matrix(sin(seq_len(60L)) + seq_len(60L) / 100, nrow = 3L,
                     dimnames = list(c("G1", "G2", "G3"), paste0("c", 1:20)))
pseudotime <- data.frame(time = seq(0, 1, length.out = 20L), row.names = colnames(expression))
paths <- c(PIDC = tempfile(fileext = ".csv"), SCODE = tempfile(fileext = ".csv"))
pidc <- run_PIDC(expression, filename = paths[["PIDC"]])
scode <- run_SCODE(expression, pseudotime, filename = paths[["SCODE"]],
                   latent_dim = 2L, n_iter = 2L, n_runs = 1L)
for (network in list(pidc, scode)) {
  stopifnot(is.data.frame(network), nrow(network) == 6L,
            identical(names(network), c("regulator", "target", "weight")),
            all(network$regulator != network$target))
}
stopifnot(all(file.exists(paths)))
unlink(paths)

check_benchmark_worker <- function() {
  previous_repo <- repo
  fixture <- tempfile("benchmark worker's ")
  dir.create(fixture)
  on.exit({
    repo <<- previous_repo
    unlink(fixture, recursive = TRUE)
  }, add = TRUE)
  repo <<- fixture
  cases <- file.path(fixture, "data/beeline/Curated/HSC", sprintf("case-%03d", 1:500))
  for (case in cases) {
    dir.create(case, recursive = TRUE)
    file.create(file.path(case, "ExpressionData.csv"))
  }
  write.csv(expression, file.path(cases[[1L]], "ExpressionData.csv"))
  write.csv(pseudotime, file.path(cases[[1L]], "PseudoTime.csv"))
  output <- file.path(fixture, "output")
  dir.create(file.path(output, "work"), recursive = TRUE)
  run_benchmark_method(output, "PPCOR", limit = 1L)
  record <- data.table::fread(file.path(output, "work/PPCOR.tsv"))
  stopifnot(nrow(record) == 1L, record$Method == "PPCOR", record$Case == "case-001",
            record$Seed == 1L, record$TimingScope == "inference_and_export",
            record$Status == "complete", record$RuntimeSeconds > 0,
            record$NetworkSHA256 == digest::digest(
              file = file.path(output, record$Network), algo = "sha256"),
            file.exists(file.path(output, "work/PPCOR_session.txt")),
            file.exists(file.path(output, "work/PPCOR_packages.tsv")),
            identical(method, "multiCSN"))
}
check_benchmark_worker()
cat("PASS R method adapters and vendor loading outside the repository working directory\n")
