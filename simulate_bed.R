library(magrittr)
library(QDNAseq)
library(data.table)

args <- commandArgs(trailingOnly = T)

# Simulates a "true" sample (a copy number profile)
.simulate.cnv <- \(n.segments, min = 1, max = 5, 
                   cna.probs = .cna.probs,
                   min.length = 30, max.length = 50, dist = min.length, genome = NULL,
                   step = NULL, left = T, normalize.sample = \(s) s) {
  if (is.null(genome)) {
    genome <- rep(2, n.segments) # the baseline is 2 in diploid organism
  }
  position <- 1
  state <- 2
  while (position < (n.segments - min.length - floor(dist / 2))) {
    state.old <- state
    if (state.old != 2) {
      while (state == state.old) {
        state <- sample(min:max, prob = cna.probs, size = 1)  
      }  
    }
    else {
      state <- sample(min:max, prob = cna.probs, size = 1) 
    }
    alteration.length <- sample(min.length:max.length, size = 1)
    end.position <- min(position + alteration.length, n.segments)
    genome[position:end.position] <- state
    position <- end.position
  }
  genome
}

compute.chrm.end <- \(n.bins) {
  bins <- getBinAnnotations(binSize = n.bins, verbose = F)
  chrm <- bins@data$chromosome
  n.chrm <- length(chrm %>% unique)
  sapply(1:(length(chrm) - 1), \(i) chrm[i] != chrm[i + 1]) %>% which
}

# finds the breakpoints from CNA data, returned as a boolean array
get.breakpoints <- \(cna.profile) {
  c(F, sapply(2:length(cna.profile), \(i) cna.profile[i] != cna.profile[i-1])) 
}

segment.df <- \(cn.profile) {
  is.breakpoint <- get.breakpoints(cn.profile)
  starts <- which(is.breakpoint)
  data.frame(
    start = c(1,starts),
    end = c(starts-1, length(cn.profile)),
    value = cn.profile[c(1,is.breakpoint %>% which)]
  )
}

.cna.probs <- c(2,2,2,2,2,2,2,2,2,2,2,1,1,3,3,3,3,3,3,4,4,5) %>% 
  {table(.) / length(.)} 

# TODO : make sure that no alterations are introduced in the region 230:300

# OBS,this only simulates in the first chromosome..
generate.bed.files <- 
  \(min.cn = 1, max.cn = 5, min.bins = 30, max.bins = 50, n.bins = 500,
    save.path = "copynumber.bed", sample.probs = .cna.probs,
    chr.nr = 1, remove = 230:300, min.segs = 3) {
    print(getwd())
    # used.bins <- readRDS(path.to.used) # used regions that we simulate alterations in
    # n.used <- length(used.bins)
    n.pos <- sum(getBinAnnotations(500, verbose = F)@data$chromosome == chr.nr) 
    ends <- c(0, compute.chrm.end(n.bins))
    oo <- options(scipen = 100, digits = 4) # set printing option for .bed
    on.exit(oo)
    # initiate the .bed files
    bed.file <- data.table(
      chrom = character(), chr_start = integer(), chrom_stop = integer(),
      num_positions = integer(), copy_number = integer()  
    )
    mult <- n.bins * 1000
    # NOTE: note now that we do not have any alterations ranging over entire chromosomes
    kept.region <- (1:n.pos)[-remove]
    n.used <- length(kept.region)
    left.inds <- (min.bins + 1):(min(remove) - 1) 
    right.inds <- (max(remove) + 1):(n.pos - min.bins - 1)
    genome <- 2 + numeric(n.pos)
    n.segs <- 0
    # make sure that we introduce enough alterations 
    while (n.segs < min.segs) {
      genome[left.inds] <- .simulate.cnv(length(left.inds))
      genome[right.inds] <-.simulate.cnv(length(right.inds))
      n.segs <- sum(segment.df(genome)$value != 2)
      if (n.segs < min.segs) {
        genome <- 2 + numeric(n.pos)
      }
    }
    df <- genome %>% segment.df
    add <- ends[chr.nr]
    for (i in 1:nrow(df)) {
      pos <- df[i, 1] - 1; end.pos <- df[i, 2]; state <- df[i, 3]
      new.row <- 
        list("1", (pos + add) * mult, (end.pos + add) * mult, (end.pos - pos + 1) * mult, state)
      bed.file <- rbind(bed.file, new.row)
    }
    format(bed.file, scientific = F)
    write.table(
      bed.file, save.path, quote = F, sep = "\t", row.names = F, col.names = T
    )
  }

save.path <- paste0("../simulated_data/bed_files2/simulation", args[1], ".bed")
generate.bed.files(save.path = save.path, chr.nr = 1)

