## ----- Compute concordance of two samples ------------------------------------

get.breakpoints <- \(cna.profile) {
  c(F, sapply(2:length(cna.profile), \(i) cna.profile[i] != cna.profile[i-1])) 
}

.compare.bp <- \(orig.bp, den.bp, margin = 2) {
  n.predictions <- length(den.bp)
  tp <- 0; fp <- 0; fn <- 0
  for (i in 1:length(orig.bp)) {
    within.margin <- abs(den.bp - orig.bp[i]) <= margin
    if (any(within.margin)) { # check is we are within the tolerance
      to.remove <- which(within.margin)[1] # the first bp
      den.bp <- den.bp[-to.remove] # remove from the set of possible breakpoints
      tp <- tp + 1
    } else {
      fn <- fn + 1
    }
  }
  fp <- n.predictions - tp
  list(tp = tp, fp = fp, fn = fn)
}

.metrics <- \(orig.bp, den.bp, margin = 2) {
  cm <- .compare.bp(orig.bp, den.bp, margin = margin)
  precision <- cm$tp / (cm$tp + cm$fp)
  recall <- cm$tp / (cm$tp + cm$fn)
  F1 <- 2 * precision * recall / (precision + recall)
  data.frame(
    precision = precision, 
    recall = recall,
    F1 = F1,
    n.true = length(orig.bp),
    n.pred = length(den.bp)
  ) 
}

## ----- Run validation for ichorCNA -------------------------------------------

run.ichorCNA.pipeline <- 
  \(path.to.pipeline = "../ichorCNA/scripts/snakemake/") {
    orig.path <- getwd()
    print(orig.path)
    setwd(path.to.pipeline)
    print(getwd())
    command <- paste0("snakemake -s ichorCNA.snakefile --cores 10")
    system(command)
    setwd(orig.path)
  }

get.ichorCNA.results <- \(
  path.to.results = "../ichorCNA/scripts/snakemake/results/ichorCNA/",
  samples = NULL, bins = 500000) {
  if (is.null(samples)) {
    samples <- list.dirs(path.to.results, full.names = F, recursive = F) 
    print(length(samples))
  }
  lapply(samples, \(sample) {
    txt.path <- paste0(path.to.results, sample, "/", sample, ".seg.txt")
    cn.results <- read.table(txt.path, header = T)
    states <- cn.results$copy.number 
    positions <- cbind(floor(cn.results[, c("start", "end")] / bins), states)
    profile <- numeric(positions$end[nrow(positions)])
    profile[1] <- positions$states[1]
    for (i in 1:nrow(positions)) {
      inds <- (positions$start[i]):(positions$end[i])
      profile[inds] <- positions$states[i]
    }
    param.path <- 
      paste0(path.to.results, sample, "/", sample, ".params.txt")
    purity <- readLines(param.path)[2] %>% strsplit("\t") %>% {.[[1]][2]} %>%
      as.numeric
    list(
      profile = profile, purity = purity, sample.name = sample, 
      start = positions$start, end = positions$end
    )
  })
}

path.to.results <- "../ichorCNA/scripts/snakemake/results/ichorCNA/"
evaluate.ichorCNA <- 
  \(ichorCNA.list = NULL, results.path = path.to.results, 
    keep = NULL, path.to.bed = "bed_files/") {
  if (is.null(ichorCNA.list)) {
    ichorCNA.list <- 
      get.ichorCNA.results(path.to.results)
  }
  lapply(ichorCNA.list, \(results) {
    bed.file <- 
      paste0(path.to.bed, strsplit(results$sample.name, "_")[[1]][1], ".bed")
    true <- get.true.values(bed.file)
    if (is.null(keep)) {
      keep <- c(1, intersect(1:length(true), results$start))
    }
    metrics <- compute.metrics(true[keep], results$profile[keep])
    df <- data.frame(
      metrics,
      method = "ichorCNA",
      purity.est = results$purity,
      sample = results$sample.name
    )
    df
  }) %>% do.call(rbind, .)
}

# TODO : evaluate on one file individually..
get.bcp.results <- \(sample.name, path.to.results = "bcp_objects/", pp = T) {
  path.to.file <- paste0(path.to.results, sample.name, "_bcp.RData")
  data <- readRDS(path.to.file) 
  if (pp) data %<>% post.process.bcp
  data
}

## ----- Run the simulation sample ---------------------------------------------

# TODO : separate this function into different parts where we can 

run.evaluation.wrapper <- 
  \(path.to.bam = "bam_files_25", path.to.bed = "bed_files") {
  bam.files <- list.files(path.to.bam, full.names = T, pattern = "\\.bam$") 
  bed.files <- list.files(path.to.bed) 
  run.ichorCNA.pipeline() # to be able to get the results
  pblapply(bam.files, \(file) {
    file.name <- strsplit(file, "/")[[1]] %>% {.[length(.)]}
    bed.file <- paste0(
      path.to.bed, "/",
      strsplit(file.name,  "[_\\.]")[[1]] %>% .[1], 
      ".bed"
    )
    run.evaluation(bam.path = file, bed.path = bed.file)
  }) %>% do.call(rbind, .)
}

##### ----- ENTROPY BASED METRICS ----------------------------------------------

compute.NMI <- \(pred, true) {
  .to.labels <- \(x) {
    segs <- segment.df(x)
    lens <- segs %$% (end - start + 1)
    rep(1:nrow(segs), times = lens)
  }
  pred.labels <- .to.labels(pred)
  true.labels <- .to.labels(true)
  NMI(pred.labels, true.labels)
}
