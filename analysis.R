library(ggplot2)
library(magrittr)
library(rstudioapi)
library(devtools)
library(QDNAseq)
library(data.table)
library(pbapply)
library(bcp)
library(yaml)
library(tidyverse)
library(pracma)
library(aricode)
library(cowplot)
library(ggpubr)
library(stringr)

# set the working directory to source file location, and print location
getActiveDocumentContext()$path %>% dirname %>% setwd %T>% print 

##  ----- Settings -------------------------------------------------------------

# plot settings
theme_set(theme_bw(base_size = 12)) # set default theme and text size
cb.palette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", 
  "#D55E00", "#CC79A7"
) # color-blind friendly palette
options(ggplot2.discrete.colour = cb.palette) # set cb.palette as default

## ----- Load required functions -----------------------------------------------

source("simulate_data.R")
source("handle_data.R") 
source("insilico_data.R")
source("evaluation.R")
source("segmentation.R")

## ----- Initial analysis : show that the method works on simulated data -------

noise.levels <- c(0.5, 1, 2, 4) # noise levels to include

dataset <- simulate.dataset(
  1, 1, seq(0.05, 0.5, by = 0.05), N.segments = 1000,
  noise.mean = c(0.5, 1, 2, 4), noise.sd = rep(0, 4),
  cna.dep.noise = T, downsample = T, min.length = 30, max.length = 50
)

# run evaluation on synthetic dataset
results.synthetic <- pblapply(1:nrow(dataset$CN), \(i) {
  den <- bcp(dataset$CN[i, ] %>% as.numeric, burnin = 500, mcmc = 1000, return.mcmc = T,
             p0 = 0.01)
  bps.pred <- bcp.pp.sd(den$data[, 2], den)
  bps.true <- which(get.breakpoints(dataset$segment[i, ] %>% as.numeric))
  .metrics(bps.true, bps.pred)
}) %>% do.call(rbind, .)

results.synthetic$purity <- dataset$purity %>% as.factor
results.synthetic$noise.level <- dataset$noise.mean %>% as.factor

# set F1 score and precision to zero if no changepoints are detected
results.synthetic[which(is.na(results.synthetic$F1)), "F1"] <- 0 
results.synthetic[which(is.na(results.synthetic$precision)), "precision"] <- 0 

# plot the results
ggplot(results.synthetic, aes(x = purity, y = F1, color = noise.level, group = purity)) + 
  geom_boxplot() + facet_wrap(noise.level ~ ., ncol = 2) +
  theme(legend.position = "top") +
  stat_summary(fun = mean, geom = "point", shape = 8, size = 2, 
               color="red", fill="red") +
  labs(color = "noise level")

##### ----- Comparing BayesCNA to ichorCNA and QDNAseq -------------------------

############################## QDNAseq #########################################

# read in the the processed data
QDNAseq <- files.to.dataset("processed_samples_new", segments = F)
dataset.QDNAseq <- QDNAseq$data %>% as.matrix

# read in QDNAseq segmentation
QDNA.seg <- files.to.dataset("processed_samples_new", segments = T)
segs <- QDNA.seg$data %>% as.matrix
QDNAseq.bps <- apply(segs, 1, \(r) which(get.breakpoints(r)))

## ----- Get the true values ---------------------------------------------------

# path to bed is the path to the folder containing generated bed files
true.profiles <- get.true.profiles(dataset.QDNAseq, path.to.bed = "bed_files2/")
true.profiles <- true.profiles[, which(QDNAseq$used)]
baseline.bps <- 
  pbapply(true.profiles, 1, \(s) which(get.breakpoints(s)), simplify = F)

## ----- Run BCP on the data processed by QDNAseq ------------------------------

denoised <- pblapply(1:nrow(dataset.QDNAseq), \(i) {
  bcp(dataset.QDNAseq[i, ] %>% as.numeric %>% log2, p0 = 0.01,
      burnin = 500, mcmc = 2000, return.mcmc = T) 
})
# compute the breakpoints identified by bcp
bcp.breakpoints <- pblapply(denoised, \(o) {
  bcp.pp.sd(o$data[,2], o, eps = 0.05) 
})

## ----- Run evaluation --------------------------------------------------------

results <- pblapply(1:ncol(dataset.QDNAseq), \(i) { 
  name <- rownames(dataset.QDNAseq)[i]
  baseline <- true.profiles[i, ]
  bps.base <- which(get.breakpoints(baseline))
  bps.bcp <- bcp.breakpoints[[i]] 
  bps.QDNAseq <- QDNAseq.bps[[i]]
  metrics.peak <- .metrics(bps.base, bps.bcp, margin = 2)
  metrics.peak$method <- "BCP"
  metrics.QDNAseq <- .metrics(bps.base, bps.QDNAseq, margin = 2)
  metrics.QDNAseq$method <- "QDNAseq"
  df.combined <- rbind(metrics.peak, metrics.QDNAseq) 
  split <- strsplit(name, "_")[[1]]
  df.combined$coverage <- split[2] %>% sapply(\(s) {
    if (s == "015") s <- "0.15x coverage"
    else if (s == "02") s <- "0.20x coverage"
    s
  })
  df.combined$purity <- split[3]
  df.combined$sample <- name
  df.combined
}) %>% do.call(rbind, .)

results$purity %<>% as.numeric

results$precision[which(is.na(results$precision))] <- 0
results$F1[which(is.na(results$F1))] <- 0

p.precision <- 
  ggplot(results, aes(x = purity, y = precision, color = method, 
                      group = interaction(method, purity))) + 
  geom_boxplot() +
  facet_grid(.~ coverage) + 
  theme(legend.position = "top") + labs(x = "") +
  scale_color_manual(values = c("#999999", "#56B4E9"))
p.recall <- 
  ggplot(results, aes(x = purity, y = recall, color = method, 
                      group = interaction(method, purity))) + 
  geom_boxplot() +
  facet_grid(.~ coverage) + 
  theme(legend.position = "top") + labs(x = "") +
  scale_color_manual(values = c("#999999", "#56B4E9"))
p.F1 <- 
  ggplot(results, aes(x = purity, y = F1, color = method, 
                      group = interaction(method, purity))) + 
  geom_boxplot() +
  facet_grid(.~ coverage) + 
  theme(legend.position = "top") + labs(x = "Purity") +
  scale_color_manual(values = c("#999999", "#56B4E9"))

ggarrange(p.precision, p.recall, p.F1, common.legend = T, ncol = 1)

############################# ICHORCNA #########################################

## run ichorCNA ----------------------------------------------------------------

# folder in which the results are saved 
ichorCNA.path <- paste0(getwd(), "/ichorCNA/scripts/snakemake/results")
ichorCNA.path <- "~/cnv-detection-local/ichorCNA/scripts/snakemake/results/"
# write a sample yaml, ichorCNA will be run on the files in the folder
write.sample.yaml(
  path.to.bam = paste0(getwd(), "/mixed_samples_new"), # bam files used in processing
  path.to.normal = NULL
)
compute.index(path.to.bam = "mixed_samples_new/") # path to produced bam files
run.ichorCNA.pipeline() # run ichorCNA
ichorCNA.results <- 
  get.ichorCNA.results(path.to.results = paste0(ichorCNA.path, "ichorCNA/")) 

## ichorCNA --------------------------------------------------------------------

# dataset processed by ichorCNA
ichorCNA.dataset <- 
  wig.to.dataset(results.path = paste0(ichorCNA.path, "/ichorCNA")) 
# remove problematic bins
to.keep <- 
  which(apply(ichorCNA.dataset, 2, \(col) sum(is.na(col)) == 0))
ichorCNA.dataset %<>% {.[, to.keep]} # remove problematic bins

## Get true values for evaluation ----------------------------------------------¨

true.profiles <- get.true.profiles(ichorCNA.dataset, path.to.bed = "bed_files2/")
true.profiles %<>% {.[, to.keep]}

baseline.bps <- 
  pbapply(true.profiles, 1, \(s) which(get.breakpoints(s)), simplify = F)

## Run BCP ---------------------------------------------------------------------

denoised <- pblapply(1:nrow(ichorCNA.dataset), \(i) {
  bcp(ichorCNA.dataset[i, ] %>% as.numeric, p0 = 0.01, 
      burnin = 500, mcmc = 2000, return.mcmc = T)
})
bcp.breakpoints <- pblapply(denoised, \(o) {
  bcp.pp.sd(o$data[, 2], o, eps = 0.05) 
})

## ichorCNA --------------------------------------------------------------------

profiles <- 
  lapply(ichorCNA.results, \(s) s$profile[to.keep]) %>% do.call(rbind, .) 
ichorCNA.bps <- apply(profiles, 1, \(r) which(get.breakpoints(r)))
ichorCNA.names <- sapply(ichorCNA.results, \(l) l$sample.name)
ichorCNA.id <- sapply(ichorCNA.names, \(name) strsplit(name, "_")[[1]][1])

## Evaluation ------------------------------------------------------------------

results <- pblapply(1:nrow(ichorCNA.dataset), \(i) { 
  name <- rownames(dataset.QDNAseq)[i]
  bps.base <- baseline.bps[[i]] 
  bps.bcp <- bcp.breakpoints[[i]]
  sample.id <- paste0(strsplit(name, "_")[[1]][1], "_0")
  bps.ichorCNA <- ichorCNA.bps[[which(ichorCNA.names == sample.id)]]
  metrics.peak <- .metrics(bps.base, bps.bcp, margin = 2)
  metrics.peak$method <- "BayesCNA"
  metrics.ichorCNA <- .metrics(bps.base, bps.ichorCNA, margin = 2) 
  metrics.ichorCNA$method <- "ichorCNA"
  df.combined <- rbind(metrics.peak, metrics.ichorCNA) 
  split <- strsplit(name, "_")[[1]]
  df.combined$coverage <- split[2] %>% sapply(\(s) {
    if (s == "015") s <- "0.15x coverage"
    else s <- "0.20x coverage"
    s
  })
  df.combined$purity <- split[3]
  df.combined$sample <- name
  df.combined
}) %>% do.call(rbind, .)

results$purity %<>% as.numeric

results$precision[which(is.na(results$precision))] <- 0
results$F1[which(is.na(results$F1))] <- 0

p.precision <- 
  ggplot(results, aes(x = purity, y = precision, color = method, 
                      group = interaction(method, purity))) + 
  geom_boxplot() +
  facet_grid(.~ coverage) +
  theme(legend.position = "top") + labs(x = "")
p.recall <- 
  ggplot(results, aes(x = purity, y = recall, color = method, 
                      group = interaction(method, purity))) + 
  geom_boxplot() +
  facet_grid(.~ coverage) + 
  theme(legend.position = "top") + labs(x = "")
p.F1 <- 
  ggplot(results, aes(x = purity, y = F1, color = method, 
                      group = interaction(method, purity))) + 
  geom_boxplot() +
  facet_grid(.~ coverage) + 
  theme(legend.position = "top") + labs(x = "Purity")

ggarrange(p.precision, p.recall, p.F1, common.legend = T, ncol = 1)

##### ----- Cell line data -----------------------------------------------------

# Available here :  https://github.com/elakatos/liquidCNA_data

# rename.files("cell-line/")
exp.data <- files.to.dataset("data/cell-line", new.format = F)
segments <- files.to.dataset("data/cell-line", segments = T, new.format = F)

## Remove unused samples -------------------------------------------------------

# remove faulty samples (no technical replicate) and F samples
exp.data %<>% {.[!(rownames(.) %in% c("E5_1", "E5_2")), ]}
segments %<>% {.[!(rownames(.) %in% c("E5_1", "E5_2")), ]}
rows.to.remove <- which(substring(rownames(exp.data), 1, 1) == "F")
exp.data %<>% {.[-rows.to.remove, ]}
segments %<>% {.[-rows.to.remove, ]}

## Remove problematic bins -----------------------------------------------------

# remove blacklisted regions by QDNAseq
exp.data %<>% remove.missing
segments %<>% remove.missing

to.blacklist <- readRDS("blacklist_500k.RData")
exp.data %<>% {.[, -to.blacklist]}
segments %<>% {.[, -to.blacklist]}

epsilon <- 0.01
has.zero <- apply(segments, 1, \(x) any(x < epsilon))
apply(segments[has.zero, ], 1, \(x) {
  inds <- which(x < epsilon)
  c(min(inds), max(inds))
}) %>% print

# remove problematic regions
remove <- c(1191, 1690:1965)
exp.data %<>% {.[, -remove]}
segments %<>% {.[, -remove]}  

## Run BCP on the experimental data --------------------------------------------

# baseline samples
baselines <- which(
  rownames(exp.data) %in% c("A1_1", "A2_2", "A3_1", "A4_1", "A5_1", "A6_1")
)

baseline.samples <- exp.data[baselines, ]
baseline.segments <- segments[baselines, ]

# TODO : only denoise the ones that we will use as ref.
baseline.bcp <- pblapply(1:nrow(baseline.samples), \(i) {
  o <- bcp(baseline.samples[i, ] %>% as.numeric, p0 = 0.01,
           burnin = 500, mcmc = 2000, return.mcmc = T) 
  bcp.pp.sd(o$data[,2], o, eps = 0.05)
})

## Get QDNAseq baseline --------------------------------------------------------

baseline.QDNAseq <- apply(baseline.segments, 1, \(x) {
  which(get.breakpoints(x))
}, simplify = F)

names(baseline.bcp) <- rownames(baseline.samples)

## Load downsampled data -------------------------------------------------------

# read the data into .GlobalEnv
load("data/Copynumber_calls_500.RData")

dataset <- {
  raw.data <- copyNumbersSegmented@assayData$copynumber %>% as.matrix %>% t
  used.bins <- copyNumbersSegmented@featureData@data$use
  raw.data[, used.bins]
}
# QDNAseq segmentation
segments.sub <- {
  raw.data <- copyNumbersSegmented@assayData$segmented %>% as.matrix %>% t
  used.bins <- copyNumbersSegmented@featureData@data$use
  raw.data[, used.bins]
}

## Remove unused samples -------------------------------------------------------

dataset %<>% {.[!(rownames(.) %in% c("E5_1-10", "E5_2-10")), ]}
segments.sub %<>% {.[!(rownames(.) %in% c("E5_1-10", "E5_2-10")), ]}
rows.to.remove <- which(substring(rownames(dataset), 1, 1) == "F")
dataset %<>% {.[-rows.to.remove, ]}
segments.sub %<>% {.[-rows.to.remove, ]}

## Remove problematic bins -----------------------------------------------------

blacklist.extended <- readRDS("blacklist_500k.RData")
dataset %<>% {.[, -blacklist.extended]}
segments.sub %<>% {.[, -blacklist.extended]}
remove <- c(1191, 1690:1965)
dataset %<>% {.[, -remove]}
segments.sub %<>% {.[, -remove]}  

## Run denoising -----------------------------------------------------

# run denoising on the deeply sequenced data
bps <- pblapply(1:nrow(dataset), \(i) { 
  chain.1 <- bcp(dataset[i, ] %>% as.numeric, p0 = 0.01,
                 burnin = 500, mcmc = 2000, return.mcmc = T) 
  bps.after <- bcp.pp.sd(chain.1$data[, 2], chain.1, eps = 0.05) 
  l <- ncol(dataset)
  genome <- numeric(l)
  bps.after <- c(1, bps.after, l)
  for (i in 1:(length(bps.after) - 1)) {
    inds <- bps.after[i]:bps.after[i + 1]
    genome[inds] <- median(chain.1$data[inds, 2], trim = 0.05)
  }
  list(breakpoints = which(get.breakpoints(genome)), profile = genome)
})

sample.bps <- lapply(bps, \(s) s$breakpoints) 
profiles <- lapply(bps, \(s) s$profile) 
profiles %<>% do.call(rbind, .)

rownames(dataset) <- sapply(rownames(dataset), \(s) strsplit(s, "-")[[1]][1])
names(sample.bps) <- rownames(exp.data)

## Get QDNAseq segmentation ----------------------------------------------------

QDNA.bps <- apply(segments.sub, 1, \(row) which(get.breakpoints(row)))
names(QDNA.bps) <- rownames(exp.data)

## Compute condordance ---------------------------------------------------------

results <- lapply(1:nrow(profiles), \(i) {
  bps <- profiles[i, ]
  sample.name <- names(sample.bps)[i]
  id <- substring(sample.name, 2, 2) %>% as.integer
  sample.id <- substring(sample.name, 1, 1)
  if (sample.id != "A") {
    id <- id + 1
  }
  baseline <- 
    construct.profile(baseline.samples[id, ] %>% as.numeric, baseline.bcp[[id]])
  df <- data.frame(NMI = NMI(baseline, bps), sample.name = sample.name)
  df$method <- "BCP"
  df
}) %>% do.call(rbind, .)

results.QDNA <- lapply(1:nrow(segments.sub), \(i) {
  bps <- segments.sub[i, ]
  sample.name <- names(QDNA.bps)[i]
  id <- substring(sample.name, 2, 2) %>% as.integer
  sample.id <- substring(sample.name, 1, 1)
  if (sample.id != "A") {
    id <- id + 1
  }
  baseline <- 
    construct.profile(segments.sub[id, ], baseline.QDNAseq[[id]])
  df <- data.frame(NMI = NMI(baseline, bps), sample.name = sample.name)
  df$method <- "QDNAseq"
  df
}) %>% do.call(rbind, .)

## Plot the results ------------------------------------------------------------

id <- c("A1", "A2", "A3", "A4", "A5", "A6", "B1", "B2", "B3", "B4", "B5",
        "C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3", "D4", "D5",
        "E1", "E2", "E3", "E4")

metadata <- data.frame(
  id = rep(id, each = 2),
  purity = c(rep(c(100, 25, 12.5, 6.25, 3.125), times = 2 * c(6, 5, 5, 5, 4)))
)
metadata <- metadata[-c(3), ] # A2 only have one replicate

results.combined <- cbind(
  rbind(results, results.QDNA), metadata
)

ggplot(results.combined, 
       aes(x = 0.81 * purity, y = NMI, group = interaction(method, purity), color = method)) + 
  geom_boxplot() + 
  theme(legend.position = "top") +
  geom_point(position = position_jitterdodge(), alpha = 0.3)  +
  xlim(c(1,23)) +
  labs(x = "Purity")
