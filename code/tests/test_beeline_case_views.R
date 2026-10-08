library(data.table)
source("code/functions/beeline_case_views.R")
for (case in c("HSC-2000-3-30", "VSC-2000-9")) {
  actual <- read_beeline_case(getwd(), case)
  expected <- if (case == "HSC-2000-3-30") c(26L, 336L, 11L) else c(15L, 531L, 11L)
  stopifnot(identical(unname(vapply(actual, nrow, integer(1))), expected),
            !anyDuplicated(actual$edges[, .(Method, regulator, target)]),
            sum(actual$summary$predicted_edges) == nrow(actual$edges),
            all(abs(actual$edges$within_network_signed_rank) <= 1))
}
# A minimal repository contains only formal inputs, no topology directory.
isolated <- tempfile("beeline-direct-")
dir.create(isolated)
files <- "results/beeline/evaluation/case_metrics.tsv"
for (case in c("HSC-2000-3-30", "VSC-2000-9")) {
  group <- sub("-.*", "", case)
  files <- c(files, file.path("data/beeline/Curated", group, case, "refNetwork.csv"),
    list.files("results/beeline/networks", paste0("__", case, "[.]csv$"),
               recursive = TRUE, full.names = TRUE))
}
for (file in files) {
  destination <- file.path(isolated, file)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file, destination))
}
stopifnot(!dir.exists(file.path(isolated, "results/beeline/topology")))
for (case in c("HSC-2000-3-30", "VSC-2000-9")) {
  stopifnot(identical(read_beeline_case(isolated, case), read_beeline_case(getwd(), case)))
}
message("PASS: both formal case views satisfy the display contract and load without topology.")
unlink(isolated, recursive = TRUE)
