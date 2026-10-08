#!/usr/bin/env Rscript
# Isolated integration tests: no formal result directories are written.
source("code/beeline/run_core_scaling.R")
repo <- normalizePath(getwd())
rscript <- file.path(R.home("bin"), "Rscript")
root <- tempfile("scaling's test ")
dir.create(root)
expect_error <- function(code) {
  failure <- tryCatch({force(code); FALSE}, error = function(e) TRUE)
  stopifnot(failure)
}

tryCatch({
  # The real worker still performs inference and excludes input/serialization.
  set.seed(27)
  input <- file.path(root, "actual input.rds")
  x <- matrix(rnorm(500*80),500,80,dimnames=list(paste0('C',1:500),paste0('G',1:80)))
  x[,2] <- 0.8*x[,1]+rnorm(500,sd=0.2)
  metadata <- data.frame(Case='actual',Scenario='cell_scaling',RequestedGenes=80,RequestedTFs=80)
  saveRDS(list(x=x,pt=seq(0,1,length.out=500),metadata=metadata),input)
  outputs <- list()
  for (cores in c(1L,2L)) {
    output <- file.path(root,paste0('actual_',cores));dir.create(output)
    row <- run_scaling_worker(input,output,cores,
      file.path(repo,'code/beeline/run_core_scaling.R'),rscript,0.01)
    stopifnot(row$Cores==cores,row$RuntimeSeconds>0,row$MemoryTreePeakMb>0,
              row$MemoryTreePeakMb>=row$MemoryTreeBeforeMb,
              row$MemorySamples>=1,row$TimingScope=='public_inference_only')
    outputs[[as.character(cores)]] <- readRDS(file.path(output,'network.rds'))
  }
  stopifnot(identical(outputs[['1']],outputs[['2']]))
  message('PASS: real R workers, RSS sampling and cross-core output identity')

  # A lightweight R worker exercises all 60 scheduler paths and publication.
  fixture <- file.path(root,'fixture');dir.create(fixture)
  for (file in c('code/beeline/run_core_scaling.R','code/beeline/summarize_scaling.R','code/functions/runtime.R')) {
    dest <- file.path(fixture,file);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    stopifnot(file.copy(file.path(repo,file),dest))
  }
  raw <- file.path(fixture,'data/beeline/scaling');dir.create(raw,recursive=TRUE)
  manifest <- data.frame(Case=paste0('case',1:10),Scenario='cell_scaling',
    RequestedCells=2,TotalGeneCount=2,RequestedGenes=2,RequestedTFs=1)
  write.csv(manifest,file.path(raw,'scaling_manifest.csv'),row.names=FALSE)
  for (case in manifest$Case) {
    d <- file.path(raw,'cell_scaling',case);dir.create(d,recursive=TRUE)
    write.csv(data.frame(C1=c(1,2),C2=c(2,3),row.names=c('G1','G2')),file.path(d,'ExpressionData.csv'))
    write.csv(data.frame(pt=c(0,1),row.names=c('C1','C2')),file.path(d,'PseudoTime.csv'))
  }
  stub <- file.path(root,'stub.R')
  writeLines(c(
    'measure_scaling_case <- function(input_file, output_dir, cores) {',
    'z <- readRDS(input_file); output <- output_dir',
    'if (file.exists(file.path(dirname(dirname(dirname(output))),"fail-worker"))) stop("injected failure")',
    'file.create(file.path(output,"ready"))',
    'while (!file.exists(file.path(output,"go"))) Sys.sleep(0.005)',
    'Sys.sleep(0.02);file.create(file.path(output,"done"))',
    'saveRDS(data.frame(regulator="G1",target="G2",weight=1),file.path(output,"network.rds"))',
    'write.csv(data.frame(Case=z$metadata$Case,Scenario="cell_scaling",AxisValue=2,Cores=cores,RuntimeSeconds=1,OutputEdges=1,TimingScope="public_inference_only"),file.path(output,"timing.csv"),row.names=FALSE)',
    '}'
  ),stub)
  fake_r <- file.path(root,'fixture Rscript')
  python <- Sys.which("python3")
  writeLines(c(paste0("#!", python), "import os,re,sys",
    paste0("expression = re.sub(r'^source\\(.*?\\);', ",
      encodeString(paste0("source(", encodeString(stub, quote='\"'), ");"), quote='\"'), ", sys.argv[-1], count=1)"),
    paste0("os.execv(",encodeString(rscript,quote='\"'),", [",encodeString(rscript,quote='\"'),", '--vanilla', '-e', expression])")), fake_r)
  Sys.chmod(fake_r,'0755')
  # Reproduce completion between the marker check and the RSS query.
  race_output <- file.path(root,'completion race');dir.create(race_output)
  real_sampler <- scaling_tree_rss
  samples <- 0L
  scaling_tree_rss <- function(handle) {
    samples <<- samples+1L
    if (samples==1L) return(real_sampler(handle))
    while (!file.exists(file.path(race_output,'done'))) Sys.sleep(0.005)
    stop('worker exited during RSS query')
  }
  row <- run_scaling_worker(input,race_output,1L,
    file.path(fixture,'code/beeline/run_core_scaling.R'),fake_r,0.005)
  scaling_tree_rss <- real_sampler
  stopifnot(row$RuntimeSeconds==1,row$MemoryTreePeakMb>0)
  message('PASS: worker completion during RSS sampling')
  target <- file.path(fixture,'results/beeline/core_scaling');dir.create(target,recursive=TRUE)
  tables <- c('core_scaling_runs.csv','edge_count_invariance.csv','one_to_ten_core_comparison.csv','scaling_summary_statistics.csv')
  for (file in tables) writeLines('previous result',file.path(target,file))
  original <- scaling_directory_hashes(target)
  backup <- file.path(root,'backups')
  run_core_scaling(fixture,target,backup,rscript=fake_r,interval=0.005)
  stopifnot(setequal(list.files(target),tables),nrow(read.csv(file.path(target,'core_scaling_runs.csv')))==60L,
    all(read.csv(file.path(target,'edge_count_invariance.csv'))$FullOutputExact),
    identical(list.files(dirname(target),all.files=TRUE,no..=TRUE),'core_scaling'))
  previous <- list.dirs(backup,recursive=FALSE,full.names=TRUE)[[1]]
  stopifnot(identical(scaling_directory_hashes(file.path(previous,'contents')),original),
            identical(jsonlite::read_json(file.path(previous,'sha256.json')),original))
  message('PASS: 60-run R scheduling, verified backup and four-table publication')

  # Input failure leaves the complete existing result untouched.
  current <- scaling_directory_hashes(target)
  manifest_path <- file.path(raw,'scaling_manifest.csv')
  unlink(manifest_path)
  expect_error(suppressWarnings(run_core_scaling(fixture,target,backup,rscript=fake_r)))
  stopifnot(identical(scaling_directory_hashes(target),current),
    identical(list.files(dirname(target),all.files=TRUE,no..=TRUE),'core_scaling'))
  write.csv(manifest,manifest_path,row.names=FALSE)

  # Sampling failure kills the started worker before its temporary output is removed.
  original_sampler <- scaling_tree_rss
  scaling_tree_rss <- function(handle) stop('injected sampling failure')
  expect_error(run_core_scaling(fixture,target,backup,rscript=fake_r))
  scaling_tree_rss <- original_sampler
  stopifnot(identical(scaling_directory_hashes(target),current),
    identical(list.files(dirname(target),all.files=TRUE,no..=TRUE),'core_scaling'))

  # Failed publication restores the prior complete target directory.
  expect_error(suppressWarnings(publish_scaling_result(file.path(root,'missing'),target,backup)))
  stopifnot(identical(scaling_directory_hashes(target),current))
  message('PASS: preparation/sampling failures preserve results; publication rollback restores target')
}, finally=unlink(root,recursive=TRUE))
message('Scaling tests passed.')
