# load the packages
library(magrittr)
library(QDNAseq) # TODO : might need some other package for processing
library(optparse)

##### ------ HELPER FUNCTION ---------------------------------------------------

# we need BAM files to use QDNAseq
generate.QDNAseq <- \(path.to.file, bin.size = 500, verbose = F) {
  bins<-getBinAnnotations(binSize = bin.size)
  binReadCounts(bins, bamfiles = path.to.file, verbose = verbose)
}

process.QDNAseq <- \(read.counts, chr = "1", verbose = T) {
  chromosomes <- c(as.character(1:22), "X", "Y")
  to.remove <- setdiff(chromosomes, chr)
  processed <- read.counts %>% 
    applyFilters(chromosomes = to.remove, verbose = verbose,
                 residual = T, blacklist = TRUE, mappability = 75) %>%
    estimateCorrection(verbose = verbose) %>% 
    correctBins(verbose = verbose) %>% 
    normalizeBins(verbose = verbose, method = "mode") %>%
    segmentBins(verbose = verbose) 
  processed
}

##### --------------------------------------------------------------------------

# Implement a function that manually normalizes the sample based on normal sample
# to avoid biases introduced by the GC-correction

## path.to.bam : the mixed samples that we should process
## path.to.normal : control sample used for normalization
process.sample <- \(path.to.bam, path.to.normal, chr = "1") {
   chromosomes <- c(as.character(1:22), "X", "Y")
   to.remove <- setdiff(chromosomes, chr)
   control <- generate.QDNAseq(path.to.normal) %>%
     applyFilters(residual = T, blacklist = TRUE, mappability = 75, 
                  chromosomes = to.remove) %>%
     estimateCorrection
   normal.fit <- as.matrix(control@assayData$fit)
   sample <- generate.QDNAseq(path.to.bam) %>%
     applyFilters(residual = T, blacklist = TRUE, mappability = 75, 
                  chromosomes = to.remove) %>%
     correctBins(fit = normal.fit) %>%
     normalizeBins(method = "mode") %>%
     segmentBins
   sample
  }

# define input and output argument
option.list <- list(
  make_option(c("-n", "--normal"), type = "character", default = NULL),
  make_option(c("-c", "--cancer"), type = "character", default = NULL),
  make_option(c("-o", "--output"), type = "character", default = "processed_samples_final")
)

opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)

path.to.cancer <- opt$cancer
path.to.normal <- opt$normal
outputdir <- opt$output

processed.sample <- 
  process.sample(path.to.bam = path.to.cancer, path.to.normal = path.to.cancer)

name <- strsplit(path.to.cancer, "/")[[1]][2]
name.split <- strsplit(name, "[.]")[[1]] # name of the sample
sample.name <- paste0(name.split[1], name.split[2], sep="")
print(sample.name)

save(processed.sample, file=paste0(outputdir, "/", sample.name, ".RData"))





