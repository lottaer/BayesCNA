# load the packages
library(magrittr)
library(QDNAseq) # TODO : might need some other package for processing
library(optparse)
library(ggplot2)
library(cowplot)

##### PROCESSING OF DATA #####

# we need BAM files to use QDNAseq
generate.QDNAseq <- \(path.to.file, bin.size = 500, verbose = F) {
  bins<-getBinAnnotations(binSize = bin.size)
  binReadCounts(bins, bamfiles = path.to.file, verbose = verbose)
}

process.QDNAseq <- \(read.counts, chr = "1", verbose = T) {
  chromosomes <- c(as.character(1:22), "X", "Y")
  to.remove <- setdiff(chromosomes, chr)
  processed <- read.counts %>% 
    applyFilters(chromosomes = to.remove, verbose = verbose) %>%
    estimateCorrection(verbose = verbose) %>% 
    correctBins(verbose = verbose) %>% 
    normalizeBins(verbose = verbose, method = "mode") %>%
    segmentBins(verbose = verbose) 
  processed
}

##### PURITY ESTIMATION #####

# plot the peaks of a density estimate
.plot.peaks <- \(density.est, modes, name = NULL) {
  data.frame(
    x = density.est$x,
    y = density.est$y
  ) %>%
    ggplot(aes(x = x, y = y)) + 
    geom_area(fill = "lightgrey") + geom_line() +
    geom_vline(xintercept = modes, linetype = "dashed", color = "firebrick") +
    labs(y = NULL, x = NULL) + 
    ggtitle(name)
}

# plot the peaks of a density estimate in a grid
.plot.peaks.grid <- \(fits.list, name) {
  pl.a <- plot_grid(plotlist = fits.list$plots, ncol = 2, align = "vh")
  title <- ggdraw() + draw_label(name)
  pl.a <- plot_grid(title, pl.a, ncol = 1,
                    rel_heights=c(0.1, 1))
  pl.a
}

# compute the modes of a density estimate
.compute.modes <- \(density, epsilon = 0.1) {
  cutoff <- density$y %>% max * epsilon
  density$x %<>% .[which(density$y>(0 + cutoff))]
  density$y %<>% .[which(density$y>(0 + cutoff))]
  density$x[which((diff(sign(diff(density$y))) == -2))] # %T>% print
}

# adjust a density s.t. the highest peak is 2
.density.adjusted <- \(denoised, adjust, lower, upper, epsilon = 0) {
  density.est <- density(denoised, from = lower, to = upper, adj = adjust)
  modes <- .compute.modes(density.est, epsilon = epsilon)
  modes.max <- which.max(density.est$y) # find index of highest peak
  max.x <- density.est$x[modes.max]
  mode.2 <- which.min(abs(modes-max.x))
  modes.y <- density.est$y[which(density.est$x %in% modes)]
  modes.adj <- modes / modes[mode.2] * 2
  density.est$x <- density.est$x / modes[mode.2] * 2
  list(density.est = density.est, modes.adj = modes.adj, modes.y = modes.y)
}

# evaluate the purity density
.eval.purity.density <- \(p, modes, modes.y) {
  if((length(modes) == 1) | (length(modes) >= 6)) {
    return(NA)
  } 
  mode.2 <- which(modes == 2)
  nb.modes <- ((2-mode.2+1):(length(modes)+2-mode.2))[-mode.2]
  modes <- modes[-mode.2]
  modes.y <- modes.y[-mode.2]
  sapply(1:(length(modes)), \(i) 
         min(((2+p*(nb.modes[i]-2))-modes[i])^2) * modes.y[i]^2) %>% 
    sum %>% sqrt
}

# plot the probabilities of purities for different adjust values
.plot.purity.probs <- \(p.fits, p, mins, adj, name) {
  p.fits$p <- p
  p.fits <- melt(p.fits,id='p')
  p.fits$adj <- rep(adj, each = length(p))
  pl.f <- ggplot(p.fits, aes(x=p, y=value, color=factor(adj)))
  if(!is.null(name)) {
    pl.f <- pl-f + geom_point(size=2, alpha=0.3, label = sprintf("%.2f"))
  } else {
    pl.f <- pl.f + geom_point(size=2, alpha=0.3)
  }
  pl.f <- pl.f + theme_bw() +
    geom_vline(xintercept = mins, linetype=c('dashed','dotted'), 
               colour=c('red','blue')) +
    labs(title = name, color = "adjust") +
    theme(legend.position = "top")
  pl.f
}

compute.purity <-
  \(denoised, p = seq(0.005, 1, length.out = 200),
    n.adj = 10, adj.min = 0.8, adj.max = 1.6, 
    adj = seq(adj.min, adj.max, length.out = n.adj), 
    show.plots = FALSE, grid.plots = 1:(min(4,n.adj)), name = NULL, 
    epsilon = 0.05, n.sd = 5, peaks.keep = 4) {
    denoised %<>% as.numeric
    # define the region to keep (otherwise too many modes introduced)
    lower <- denoised %>% {mean(.) - n.sd * sd(.)}
    upper <- denoised %>% {mean(.) + n.sd * sd(.)}
    fits.list <- lapply(adj , \(a) {
      d.adj <- .density.adjusted(denoised, a, lower, upper, epsilon = epsilon)
      pl <- .plot.peaks(d.adj$density.est, d.adj$modes.adj) + 
        ggtitle(sprintf("Adjust = %.1f", a)) + 
        theme(legend.position = "none", plot.title = element_text(size = 10))
      p.fits <- sapply(p, \(p) .eval.purity.density(p, d.adj$modes.adj, 
                                                    d.adj$modes.y))
      list(p.fits = p.fits, plots = pl)
    }) %>% purrr::transpose(.)
    p.fits <- fits.list$p.fits %>% as.data.frame
    p.cols <- p.fits[ ,1:length(adj)] %>% as.matrix
    which.mins <- apply(p.cols, 2, which.min) %>% unlist
    mins <- c(mean(p[which.mins], na.rm = TRUE),
              median(p[which.mins], na.rm = TRUE))
    names(mins) <- c("mean", "median")
    if(show.plots){
      .plot.peaks.grid(fits.list[grid.plots], name) %>% print
      .plot.purity.probs(p.fits, p, mins, adj, name) %>% print
    }
    if(length(mins) == 0) {
      mins <- c(NA, NA)
    }
    return(mins %>% as.data.frame)
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

get.breakpoints <- \(cna.profile) {
  c(F, sapply(2:length(cna.profile), \(i) cna.profile[i] != cna.profile[i-1])) 
}

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

# define input and output argument
option.list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL),
  make_option(c("-b", "--bed"), type = "character", default = NULL)
)

opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)

input.bam <- opt$input
input.bed <- opt$bed

# load the simulated data
qdna.obj <- generate.QDNAseq(input.bam) %>% process.QDNAseq(verbose = F) 
use <- which(qdna.obj@featureData@data$use)
processed.data <- qdna.obj@assayData$copynumber[use] * 2

true.profile <- get.true.values(input.bed)[use] 
true.segment.df <- true.profile %>% segment.df
n.peaks <- unique(true.profile) %>% length # true number of peaks

data <- numeric(nrow(true.segment.df))
# loop over the true segments and extract the data
for (i in 1:nrow(true.segment.df)) {
  inds <- true.segment.df[i, ] %$% start:end
  data[i] <- median(processed.data[inds]) # get plausible value for the segment level
}

# adjust a density s.t. the highest peak is 2
.density.adjusted <- \(denoised, adjust, lower, upper, true.values = true.profile, epsilon = 0) {
  density.est <- density(denoised, from = lower, to = upper, adj = adjust)
  density.true <- density(true.values)
  # find the number that correspons to the CN = 2 peak by finding the one with the smallest error
  CN.2.peak <- which.min(abs(density.true$y - 2))
  modes <- .compute.modes(density.est, epsilon = epsilon)
  modes.max <- which.max(density.est$y) # find index of highest peak
  max.x <- density.est$x[modes.max]
  mode.2 <- which.min(abs(modes-max.x))
  modes.y <- density.est$y[which(density.est$x %in% modes)]
  modes.adj <- modes / modes[mode.2] * 2
  density.est$x <- density.est$x / modes[mode.2] * 2
  list(density.est = density.est, modes.adj = modes.adj, modes.y = modes.y)
}

# TODO : make sure that we place the peaks in the correct positions

# TODO : adjust the peaks so that the distribution is correct...
test.adjust <- seq(0.05, 10, length.out = 1000)
peaks <- sapply(test.adjust, \(adj) density(data, adj = adj) %>% 
                  .compute.modes(epsilon = 0.01) %>% length %>% {. == n.peaks})
min.adj <- test.adjust[which(peaks)[1]]
max.adj <- which(peaks) %>% {test.adjust[.[length(.)]]}

purity <- compute.purity(data, epsilon = 0.01, 
                         adj.min = min.adj, adj.max = max.adj, show.plots = F)[2, ]

cat(as.numeric(purity)); invisible(purity)



