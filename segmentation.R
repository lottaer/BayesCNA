##### ------ Helper functions --------------------------------------------------

construct.profile <- \(sample, bps, genome.length = ncol(exp.data)) {
  genome <- numeric(genome.length)
  bps <- c(1, bps, genome.length) 
  for (i in 1:(length(bps) - 1)) {
    inds <- bps[i]:bps[i + 1]
    genome[inds] <- median(sample[inds], trim = 0.025)
  }
  genome
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

# finds the breakpoints from CNA data, returned as a boolean array
get.breakpoints <- \(cna.profile) {
  c(F, sapply(2:length(cna.profile), \(i) cna.profile[i] != cna.profile[i-1])) 
}

# construct a segment data frame from set breakpoints
seg.from.bps <- \(cn.sample, breakpoints) {
  values <- sapply(1:nrow(breakpoints), \(i) {
    cn.sample[breakpoints$start[i]:breakpoints$end[i]] %>% median
  })
  data.frame(
    start = breakpoints$start,
    end = breakpoints$end,
    value = values
  )
}

df.to.sequence <- \(segment.df) {
  lens <- segment.df %$% {end - start + 1}
  rep(segment.df$value, times = lens)
}

get.profile <- \(segment.df) {
  segment.df %$% rep(value, times = end - start + 1)
}

##### ------ Segmentation ------------------------------------------------------

find.peaks <- \(U, eps = 0.1, min.bins = 10) {
  probs <- colSums(U) / nrow(U)
  probs.normalized <- (probs - min(probs)) / (max(probs) - min(probs)) 
  probs.normalized[probs.normalized < eps] <- 0
  findpeaks(probs.normalized, minpeakdistance = min.bins, nups = 1)[, 2]
}

remove.small.segments <- \(cn.profile, n.sd = 0.5) {
  bps <- get.breakpoints(cn.profile)
  values <- cn.profile[bps] # returns the segment values
  indices <- which(bps)
  # do merge
  seg.sd <- sd(values)
  new.values <- values
  for (i in seq_along(values)[-1]) {
    if (abs(values[i] - values[i - 1]) < n.sd * seg.sd) {
      new.values[i] <- new.values[i - 1]
    }
  }
  # Reconstruct the full function to match input length
  reconstructed.func <- cn.profile
  for (i in seq_along(indices)) {
    start <- indices[i]
    end <- if (i < length(indices)) indices[i + 1] - 1 else length(cn.profile)
    reconstructed.func[start:end] <- new.values[i]
  }
  reconstructed.func
}

bcp.pp.sd <- \(sample, bcp.results, eps = 0.05, n.sd = 0.05) {
  burnin <- bcp.results$burnin
  U <- bcp.results$mcmc.rhos %>% t
  U %<>% {.[(burnin + 1):nrow(.), ]} 
  peaks <- sort(c(1, find.peaks(U, eps = eps) + 1, length(sample)))
  pp.profile <- numeric(length(sample))
  for (i in 1:(length(peaks) - 1)) {
    pp.profile[peaks[i]:peaks[i + 1]] <- median(sample[peaks[i]:peaks[i + 1]])
  }
  pp.profile.cleaned <- remove.small.segments(pp.profile, n.sd = n.sd)
  which(get.breakpoints(pp.profile.cleaned))
}

##### ------ BayesCNA ----------------------------------------------------------

run.BayesCNA <- 
  \(data, p = 0.01, eps = 0.05, eta = 0.05, n.mcmc = 2000, n.burnin = 500) {
  denoised <- bcp(data, p0 = p, burnin = n.burnin, mcmc = n.mcmc, return.mcmc = T)
  postprocessed <- bcp.pp.sd(denoised$data[, 2], denoised, eps = eps, n.sd = eta)
  profile <- construct.profile(data, postprocessed, genome.length = length(data))
  list(partition = postprocessed, profile = profile)
}

compute.CI <- \(denoised) {
  burnin <- denoised$burnin # number of burnin iterations
  iters <- denoised$mcmc # number of mcmc iterations
  means <- denoised$mcmc.means[1:iters + burnin, ] %>% t # remove burnin iterations
  summary <- apply(means, 1, \(r) {
    q <- quantile(r, probs = c(0.025, 0.975))
    data.frame(
      ymin = q[1], 
      ymax = q[2],
      cn = mean(r)
    )
  }, simplify = F) %>% do.call(rbind, .)
  summary$bin <- 1:nrow(means)
  summary
}

# plot the posterior mean and 95% CI of the mean
plot.CI <- \(denoised) {
  summary <- compute.CI(denoised)
  data <- denoised$data[, 2] # the input data
  CNA.df <- data.frame(bin = 1:length(data), cn = data)
  ggplot(summary, aes(x = bin, y = cn)) + 
    geom_point(data = CNA.df, aes(x = bin, y = cn), color = "lightgray") +
    geom_line(color = "firebrick") + 
    geom_ribbon(data = summary, mapping = aes(ymin = ymin, ymax = ymax), fill = "firebrick", alpha = 0.4) 
}

