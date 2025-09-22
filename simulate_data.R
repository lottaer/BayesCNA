library(ggplot2)
library(magrittr)

# plot settings
theme_set(theme_bw(base_size = 14)) # set default theme and text size
cb.palette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", 
  "#D55E00", "#CC79A7"
) # color-blind friendly palette
options(ggplot2.discrete.colour = cb.palette) # set cb.palette as default

.cna.probs <- c(0.09090909, 0.50000000, 0.27272727, 0.09090909, 0.04545455) 

# ----- Helper functions -------------------------------------------------------

# Simulates the baseline
.simulate.cnv <- \(n.segments, min = 1, max = 5, 
                   cna.probs = .cna.probs,
                   min.length = 12, max.length = 130, genome = NULL,
                   step = NULL, normalize.sample = \(s) s) {
  if (is.null(genome)) {
    genome <- rep(2, n.segments) # the baseline is 2 in diploid organism
  }
  position <- 1
  while (position < n.segments) {
    state <- sample(min:max, prob = cna.probs, size = 1)
    alteration.length <- sample(min.length:max.length, size = 1)
    end.position <- min(position + alteration.length, n.segments)
    genome[position:end.position] <- state
    position <- end.position
  }
  normalize.sample(genome)
}

# Simulate segment data and add random noise 
# N.samples: number of samples to simulate of each provided purity
# N.segments: number of genomic bins
# min: minimim copy number
# max: maximum copy number
# cancer.prop: purity of sample
# noise.mean: mean of Gaussian noise
# noise.sd: standard deviation of Gaussian noise
# simulated.data: data to alter
# cna.dep.noise: makes noise dependent on copy numbers
# min.length: minimal length of alterations
# max.length: maximal length of alterations
# downsample: simulates binning 500kb bins into 50kb bins
.simulate.noisy <- \(N.samples, N.segments, min = 1, max = 5,
                     cancer.prop = 1, noise.mean = 0.2, 
                     noise.sd = 0.005, simulated.data = NULL,
                     cna.dep.noise = FALSE,
                     min.length = 12, max.length = 130, downsample = FALSE) {
  if (is.null(simulated.data)) {
    simulated.data <- .simulate.cnv(N.segments, max = max, min = min,
                                    min.length = min.length, 
                                    max.length = max.length) 
  }
  # reduce the signal (normal contamination)
  segment.data <- sapply(1:N.samples, \(i) simulated.data) * cancer.prop  + 
    (1 - cancer.prop) * rep(2, N.samples) 
  if(downsample) {
    segment.data.down <- segment.data
    segment.data %<>% apply(2, rep, each = 10)
  }
  noisy <- matrix(0, nrow(segment.data), ncol(segment.data))
  if (cna.dep.noise) {
    for (i in 1:N.samples) {
      sigma0 <- max(0.01, rnorm(1, noise.mean, noise.sd)) # use a minimal noise
      sigma <- 0.2*sigma0 + 0.015*sigma0*segment.data[,i] # TODO : doible check this!!!
      noisy[, i] <- segment.data[, i] + rnorm(nrow(segment.data), 0, sigma)
    }
  }
  else {
    # create noisy measurements
    for (i in 1:N.samples) {
      sigma <- max(0.05, rnorm(1, noise.mean, noise.sd)) # use a minimal noise
      noisy[, i] <- segment.data[, i] + rnorm(nrow(segment.data), 0, sigma)
    }  
  }
  if(downsample) {
    noisy <- sapply(1:N.samples, \(s) sapply(1:N.segments, \(i) 
                                             mean(noisy[(10*i-9):(10*i),s])))
    segment.data <- segment.data.down
  }
  list(segment = segment.data, CN = noisy)
}

# Helper function, simulates dataset of one purity 
.simulate.noisy.dataset <- \(N.samples, N.segments, N.replicates = 2, 
                             max = 5, min = 1, cancer.prop = 1,
                             noise.mean = 0.2, noise.sd = 0.005,
                             simulated.data = NULL, cna.dep.noise = FALSE,
                             min.length = 12, max.length = 130, 
                             downsample = FALSE) {
  lapply(1:N.samples, \(i) .simulate.noisy(N.samples = N.replicates, 
                                           N.segments = N.segments, 
                                           max = max, min = min, 
                                           cancer.prop = cancer.prop, 
                                           noise.mean = noise.mean, 
                                           noise.sd = noise.sd, 
                                           simulated.data = simulated.data,
                                           cna.dep.noise = cna.dep.noise,
                                           min.length = min.length, 
                                           max.length = max.length, 
                                           downsample = downsample)) %>%
    {list(
      CN = lapply(., \(list) data.frame(list$CN) %>% t) %>% 
        do.call(rbind, .) %>% t,
      segment = lapply(., \(list) data.frame(list$segment) %>% t) %>% 
        do.call(rbind, .) %>% t
    )}
}

# ----- Simulate synthetic data -----------------------------------------------

# Simulates dataset of multiple purities and different noise levels
# N.samples: number of samples to simulate of each provided purity an noise level
# N.segments: number of genomic bins
# min: minimim copy number
# max: maximum copy number
# cancer.prop: vector of purities
# noise.mean: vector of means of Gaussian noise
# noise.sd: vector of standard deviation of Gaussian noise
# simulated.data: data to alter
# cna.dep.noise: makes noise dependent on copy numbers
# min.length: minimal length of alterations
# max.length: maximal length of alterations
# downsample: simulates binning 500kb bins into 50kb bins
simulate.dataset <- \(N.each, N.replicates, cancer.prop,
                      N.segments = 1000, min.CN = 1, max.CN = 5,
                      simulated.data = NULL, noise.mean = 0.2,
                      noise.sd = 0.005, cna.dep.noise = FALSE,
                      min.length = 12, max.length = 130, downsample = FALSE) {
  # the dataset will be of size N.each * length(props)
  if (length(cancer.prop) == 1) {
    cancer.prop %<>% c # add to vector to be able to apply over it
  }
  if (length(noise.mean) != length(noise.sd)) {
    stop("noise.mean and noise.sd must be of equal length")
  }
  if (length(noise.mean) == 1) {
    noise.mean %<>% c
    noise.sd %<>% c
  }
  if (length(noise.mean) != length(noise.sd)) {
    stop()
  }
  lapply(1:length(noise.mean), \(i) {
    lapply(cancer.prop,
           \(prop) {
             dataset <- .simulate.noisy.dataset(N.each, N.segments, N.replicates,
                                                cancer.prop = prop,
                                                noise.mean = noise.mean[i],
                                                noise.sd = noise.sd[i],
                                                cna.dep.noise = cna.dep.noise,
                                                min.length = min.length,
                                                max.length = max.length, 
                                                min = min.CN, max = max.CN, 
                                                downsample = downsample)
             list(CN = dataset$CN, segment = dataset$segment)
           }
    ) %>% purrr::transpose(.) %>%
      {list(
        CN = .$CN %>% do.call(cbind, .) %>% t,
        segment = .$segment %>% do.call(cbind, .) %>% t,
        purity = rep(cancer.prop, each = N.each * N.replicates),
        noise.mean = rep(
          noise.mean[i], each = N.each * N.replicates * length(cancer.prop)),
        noise.sd = rep(
          noise.sd[i], each = N.each * N.replicates * length(cancer.prop))
      )}
  }) %>% purrr::transpose(.) %>%
    {list(
      CN = .$CN %>% do.call(rbind, .),
      segment = .$segment %>% do.call(rbind, .),
      purity = .$purity %>% unlist, # concatenate the purities
      noise.mean = .$noise.mean %>% unlist,
      noise.sd = .$noise.sd %>% unlist
    )}
}


