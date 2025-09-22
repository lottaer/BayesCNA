## ----- Helper functions ------------------------------------------------------

# computes how much more weight there is in the cancer sample
compute.weight.factor <- \(cna.profile) {
  cancer <- sum(cna.profile) # this is the ground truth
  healthy <- 2 * length(cna.profile)
  cancer / healthy
}

# mix together the healthy blood data and the (simulated?) CNV data?
mix.samples <- \(X1, X2, purity) {
  alpha <- compute.weight.factor(X1)
  alpha / purity * X1 + (1 - purity) * X2
}

## ----- QDNAseq functions -----------------------------------------------------

# we need BAM files to use QDNAseq
generate.QDNAseq <- \(path.to.file, bin.size = 500, verbose = F) {
  bins<-getBinAnnotations(binSize = bin.size)
  binReadCounts(bins, bamfiles = path.to.file, verbose = verbose)
}

process.QDNAseq <- \(read.counts, chr = "1", remove.chr = T, verbose = T) {
  chromosomes <- c(as.character(1:22), "X", "Y")
  if (remove.chr) {
    to.remove <- setdiff(chromosomes, chr)
    processed <- read.counts %>% 
      applyFilters(chromosomes = to.remove, verbose = verbose,
                   residual = T, blacklist = TRUE, mappability = 75) %>%
      estimateCorrection(verbose = verbose) %>% 
      correctBins(verbose = verbose) %>% 
      normalizeBins(verbose = verbose, method = "mode") %>%
      segmentBins(verbose = verbose) 
  }
  else {
    processed <- read.counts %>% 
      applyFilters(verbose = verbose,
                   residual = T, blacklist = TRUE, mappability = 75) %>%
      estimateCorrection(verbose = verbose) %>% 
      correctBins(verbose = verbose) %>% 
      normalizeBins(verbose = verbose, method = "mode") %>%
      segmentBins(verbose = verbose) 
  }
  processed
}

get.used.bins <- \(obj, chr = "1") {
  meta.data <- obj$meta.data
  to.keep <- meta.data$chromosome == chr
  meta.data$use[to.keep]
}

remove.blacklist <- \(obj, bin.size = 500, chr = "1") {
  meta.data <- obj$meta.data
  to.keep <- meta.data$chromosome == chr
  meta.data %<>% {.[to.keep, ]}
  use <- meta.data$use
  obj$copy.number <- obj$cn.raw[use]
  obj$segmented <- obj$seg.raw[use]
  obj
}

# removes blacklisted regions from the true copy number profile
remove.unused.regions <- \(vec, obj) {
  to.keep <- obj$meta.data$use 
  vec[to.keep] %>% {.[!is.na(.)]}
}

# get the used bins
get.data <- \(qdna.obj) {
  qdna.obj@assayData$copynumber[qdna.obj@featureData@data$use]
}

## ----- Generate .BED-files to use for simulation -----------------------------

compute.chrm.end <- \(n.bins, verbose = verbose) {
  bins <- getBinAnnotations(binSize = n.bins, verbose = verbose)
  chrm <- bins@data$chromosome
  n.chrm <- length(chrm %>% unique)
  sapply(1:(length(chrm) - 1), \(i) chrm[i] != chrm[i + 1]) %>% which
}

# write wrapper around this?
.cna.probs <- c(2,2,2,2,2,2,2,2,2,2,2,1,1,3,3,3,3,3,3,4,4,5) %>% 
  {table(.) / length(.)} 
generate.bed.files <- 
  \(min.cn = 1, max.cn = 5, min.bins = 40, max.bins = 60, n.bins = 500,
    save.path = "copynumber.bed", sample.probs = .cna.probs, chr.nr = NULL) {
  oo <- options(scipen = 100, digits = 4) # set printing option for .bed
  on.exit(oo)
  ends <- c(0, compute.chrm.end(n.bins))
  # initiate empty bed file
  bed.file <- data.table(
    chrom = character(),
    chr_start = integer(),
    chrom_stop = integer(),
    num_positions = integer(),
    copy_number = integer()  
  )
  mult <- n.bins * 1000
  # NOTE: note now that we do not have any alterations ranging over entire chromosomes
  full.genome <- is.null(chr.nr)
  n.chr <- ifelse(full.genome, length(ends) - 1, 1)
  for (chr in 1:n.chr) {
    pos <- 1
    n.pos <- 
      ifelse(full.genome, ends[chr + 1] - ends[chr] + 1, ends[chr.nr + 1] - ends[chr.nr] + 1)
    while (pos < n.pos) {
      state <- sample(min.cn:max.cn, size = 1, prob = sample.probs)
      cnv.length <- sample(min.bins:max.bins, size = 1)
      end.pos <- min(pos + cnv.length, n.pos)
      if (state != 0) { 
        chr.str <- ifelse(full.genome, paste0("chr", chr), paste0("chr", chr.nr))
        add <- ifelse(full.genome, ends[chr], ends[chr.nr])
        new.row <- 
          list(chr.str, (pos + add) * mult, (end.pos + add) * mult, (end.pos - pos + 1) * mult, state)
        bed.file <- rbind(bed.file, new.row)
      }
      pos <- end.pos
    }
  }
  format(bed.file, scientific = F)
  write.table(
    bed.file, save.path, quote = F, sep = "\t", row.names = F, col.names = T
  )
}

# the path should be directed to a BAM file on disk 
get.true.values <- \(path.to.file, qdna.obj, bin.size = 500, chr = 1,
                     verbose = F) {
  bins <- getBinAnnotations(binSize = bin.size, verbose = verbose)
  cn.table <- read.table(path.to.file, header = T)
  end <- bins@data$end
  chromosomes <- bins@data$chromosome
  chr1 <- which(chromosomes == chr)
  # add <- length(chr1)
  add <- 0
  genome.binned <- numeric(length(end)) 
  for (i in 1:nrow(cn.table)) {
    curr <- cn.table[i, ]
    curr.start <- {curr$chr_start / (bin.size * 1000)} %>% floor
    curr.end <- {curr$chrom_stop / (bin.size * 1000)} %>% ceiling
    genome.binned[(add + curr.start):(add + curr.end)] <- curr$copy_number
  }
  if (!is.null(chr)) {
    chromosomes <- bins@data$chromosome
    genome.binned %<>% {.[chromosomes == chr]}
  }
  genome.binned
}


