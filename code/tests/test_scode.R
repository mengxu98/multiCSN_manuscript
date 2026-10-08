source("code/functions/runtime.R")
root <- tempfile("SCODE worker's ")
dir.create(root)
tryCatch({
  set.seed(10L)
  expression <- matrix(rnorm(5L * 30L), 5L, 30L)
  expression_path <- file.path(root, "expression.tsv")
  time_path <- file.path(root, "time.tsv")
  write.table(expression, expression_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(data.frame(cell = 0:29, time = seq(0, 1, length.out = 30)), time_path,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  vendor <- normalizePath("code/vendor/SCODE/SCODE.R")
  # Pin the seed in the test adapter only; production SCODE RNG policy is unchanged.
  adapter <- file.path(root, "seeded vendor.R")
  writeLines(c(paste0("source(", encodeString(vendor, quote = '"'), ")"), "set.seed(42L)"), adapter)
  arguments <- list(fdata = expression_path, ftime = time_path, dir = file.path(root, "actual"),
                    tfnum = 5L, pnum = 2L, cnum = 30L, maxite = 5L)
  processx::run(file.path(R.home("bin"), "Rscript"),
    c("-e", r_function_expression(adapter, "run_scode", arguments)))
  source(adapter)
  arguments$dir <- file.path(root, "reference")
  do.call(run_scode, arguments)
  expected <- sort(list.files(arguments$dir))
  stopifnot(identical(expected, sort(list.files(file.path(root, "actual")))))
  for (name in expected) stopifnot(identical(readLines(file.path(arguments$dir, name)),
                                           readLines(file.path(root, "actual", name))))
  cat("PASS SCODE isolated function call and seeded equivalence to the original vendor function\n")
}, finally = unlink(root, recursive = TRUE))
