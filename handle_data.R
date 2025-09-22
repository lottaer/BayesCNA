library(QDNAseq)

# ----- Helper functions -------------------------------------------------------

rename.files <- \(path.to.files, copy.numbers = T, ending = ".RData", pattern = ending) {
  if (ending == ".RData") {
    file.names <- list.files(path.to.files) %>% get.file.names(copy.numbers)  
  }
  else {
    file.names <- list.files(path.to.files, pattern = ending)
  }
  print(file.names)
  for (file in file.names) {
    split <- strsplit(file, "[.]")[[1]]
    print(split)
    # if we have a . in the name, we will rename the file
    if (length(split) > 2) {
      print(paste0("Re-naming: ", file))
      if (ending == ".RData") {
        new.name <- paste0(path.to.files, "/", split[1], "_", split[2], 
                           collapse = "", ending)  
      }
      else {
        new.name <- paste0(path.to.files, "/", split[1], "", split[2], 
                           collapse = "", ending)  
      }
      file.rename(paste0(path.to.files, "/", file), new.name)
    }
  }
}

# Get the name of either all copy-number files, or readcount files
get.file.names <- \(files, copy.numbers = T) {
  file.start <- ifelse(copy.numbers, "Copynumber", "Readcounts")
  to.keep <- sapply(files, \(file) strsplit(file, "_")[[1]][1] == file.start)
  files[to.keep]
}

# Removes NA:s from a dataset
remove.missing <- \(df) {
  if(is.null(dim(df))){
    df %<>% as.matrix %>% t
  } else if(dim(df)[1] == 1) {
    df %<>% as.matrix
  }
  to.remove <- apply(df, 2, \(col) any(is.na(col)))
  df[, !to.remove]
}

# ----- Functions for handling QDNAseq-objects ---------------------------------

# Gets the raw count matrix from a QDNAseq object
get.counts <- \(read.counts) {
  read.counts$readCountsFiltered@assayData$counts 
}

# Gets the CN-matrix from a QDNAseq object
get.CNVs <- \(copy.numbers) {
  copy.numbers$copyNumbersSegmented@assayData %>% 
    {list(CN = .$copynumber, segment = .$segmented)}
}

# ----- Files for loading QDNAseq data into .GlobalEnv -------------------------

load.data <- \(path.to.files = "data/cell-line") {
  .data.files <- list.files(path = path.to.files, pattern = "\\.RData$",
                            full.names = TRUE)
  .N.loaded <- 0
  for (.file in .data.files) {
    .file.name <- tools::file_path_sans_ext(basename(.file))  
    if (!exists(.file.name)) {
      .N.loaded <- .N.loaded + 1
      assign(.file.name, load(.file, .GlobalEnv) %>% mget(., parent.frame()),
             envir = .GlobalEnv)
    }
    else {
      print(paste0(.file, " already loaded. Continuing..."))
    }
  }
  all.variables <- ls(envir = .GlobalEnv)
  if ("copyNumbersSegmented" %in% all.variables) {
    rm("copyNumbersSegmented", envir = .GlobalEnv)
  }
  if ("readCountsFiltered" %in% all.variables) {
    rm("readCountsFiltered", envir = .GlobalEnv)
  }
}

files.to.dataset <- 
  \(path.to.files, copy.numbers = T, segments = F, all.file.names = NULL,
    new.format = T) {
    if (is.null(all.file.names)) {
        file.names <- list.files(path.to.files)
    }
    if (!new.format) file.names %<>% get.file.names(copy.numbers)
    object.names <- sapply(file.names, \(file) strsplit(file, "[.]")[[1]][1])
    loaded.objects <- ls(envir = .GlobalEnv)
    is.loaded <- sapply(object.names, \(obj) obj %in% loaded.objects)
    if (!all(is.loaded)) {
      print("Files not in .GlobalEnv. Loading files...")
      load.data(path.to.files)
    }
    if (new.format) {
      tmp <- get(object.names[1])
      used <- tmp$processed.sample@featureData@data$use
    }
    if (copy.numbers) {
      dataset <- object.names %>%
        pblapply(\(file) {
          if (new.format) {
            data <- get(file)$processed.sample
            tmp <- data@assayData
            use  <- data@featureData@data$use
            if (segments) {
              tmp$segmented[use] %>% t
              # tmp$segmented %>% t
            }
            else {
              tmp$copynumber[use] %>% t
              # tmp$copynumber %>% t
            }
          }
          else {
            data <- get.CNVs(get(file)) #$CN %>% t
            if (segments) {
              data$segment %>% t
            }
            else  {
              data$CN %>% t
            }
          }
        }) %>% 
        do.call("rbind", .) %>% data.frame   
    } else {
      dataset <- object.names %>%
        pblapply(\(file) get.counts(get(file)) %>% t) %>% 
        do.call("rbind", .) %>% data.frame   
    }
    # remove the files that we did not have prior to running the function
    if (!all(is.loaded)) {
      print("Cleaning loaded files...")
    }
    to.remove <- setdiff(ls(envir = .GlobalEnv), loaded.objects)
    rm(list = to.remove, envir = .GlobalEnv)
    if (new.format) {
      row.names(dataset) <- sapply(file.names, \(f) strsplit(f, "[.]")[[1]][1])
      return(list(data = dataset, used = used))
    }
    dataset
  }

add.metadata.columns <- \(df.named) {
  sample.names <- row.names(df.named) %>% sapply(\(f) strsplit(f, "[.]")[[1]][1])
  sample.names.split <- lapply(sample.names, \(s) strsplit(s, "_")[[1]][2:3]) 
  cols <- sapply(sample.names.split, \(metadata) {
    coverage.str <- metadata[1]
    ch.array <- str_split_1(coverage.str, "")
    coverage <- paste(c(ch.array[1], ".", ch.array[-1]), collapse = "")
    as.numeric(c(coverage, metadata[2])) 
  }) %>% t %>% as.data.frame
  colnames(cols) <- c("coverage", "purity")
  cbind(df.named, cols)
}

get.samples <- \(df, coverage = NULL, purity = NULL) {
  keep.covs <- keep.purs <- 1:nrow(df)
  if (!is.null(coverage)) {
    covs <- df$coverage
    keep.covs <- which(covs %in% coverage)
  }
  if (!is.null(purity)) {
    purs <- df$purity
    keep.purs <- which(purs %in% purity)
  }
  to.keep <- intersect(keep.covs, keep.purs)
  df[to.keep, ]
}

remove.metadata <- \(df) {
  df[, -which(colnames(df) %in% c("coverage", "purity"))]
}

get.true.profiles <- \(df.named, path.to.bed = "bed_files/", to.remove = NULL) {
  # TODO : extract the id's
  sample.names <- rownames(df.named)
  true <- lapply(sample.names, \(s) {
    bed.file <- paste0(strsplit(s, "_")[[1]][1], ".bed")
    get.true.values(paste0(path.to.bed, bed.file))
  }) %>% do.call(rbind, .)
  if (!is.null(to.remove)) {
    true <- true[, -to.remove]
  }
  true
}

# ----- Load ichorCNA processed data -------------------------------------------

wig.to.dataset <- \(results.path = "../ichorCNA/scripts/snakemake/results/ichorCNA") {
  files <- 
    list.files(results.path, pattern = "\\correctedDepth.txt$", recursive = T, full.names = T)
  dataset <- lapply(files, read.wig.file) %>% do.call(rbind, .)
  rownames(dataset) <- strsplit(files, "/") %>% sapply(\(l) l[length(l)]) %>%
    strsplit("[.]") %>% sapply(\(l) l[1])
  dataset
}

read.wig.file <- \(path.to.wig = "~/cnv-detection-local/ichorCNA/scripts/snakemake/results/ichorCNA/simulation1006_0/simulation1006_0.correctedDepth.txt") { 
  table <- read.table(path.to.wig, header = T)
  start <- floor(table$start / 500000)
  end <- table$end / 500000
  depth.corr <- table$log2_TNratio_corrected
  seq <- rep(NA, max(end))
  seq[start] <- depth.corr
  seq
}

# TODO : check if NA's are always introduced in the same positions
get.ichorCNA.bins <- \(path.to.wig ="../ichorCNA/scripts/snakemake/results/eval/simulation1006_0/simulation1006_0.correctedDepth.txt", bins = 500 * 1e3) {
  table <- read.table(path.to.wig, header = T)
  floor(table$start / bins)
}

write.sample.yaml <- 
  \(path.to.bam = "/Users/lotta/cnv-detection-local/cnv-detection/ichorCNA-test", 
    path.to.normal = "/Users/lotta/cnv-detection-local/cnv-detection/ichorCNA-normal", 
    bam.files = NULL,
    save.path = "../ichorCNA/scripts/snakemake/config/samples.yaml") {
    if (is.null(bam.files)) {
      bam.files <- list.files(path.to.bam, full.names = T, pattern = "\\.bam$") 
    }
    else {
      bam.files <- sapply(bam.files, \(f) paste0(path.to.bam, "/", f))
    }
    file.names <- sapply(bam.files, \(file) {
      full.name <- strsplit(file, "/")[[1]] %>% {.[length(.)]} 
      strsplit(full.name, "[.]")[[1]][1]
    })
    if (is.null(path.to.normal)) {
      yaml.list <- setNames(as.list(bam.files), file.names)
      yaml.structure <- list(samples = yaml.list)
      write_yaml(yaml.structure, save.path)  
    }
    else {
      normal.files <- list.files(path.to.normal, full.names = T, pattern = "\\.bam$")
      normal.file.names <- sapply(normal.files, \(file) {
        full.name <- strsplit(file, "/")[[1]] %>% {.[length(.)]} 
        strsplit(full.name, "[.]")[[1]][1]
      })
      all.files <- c(bam.files, normal.files)
      all.file.names <- c(file.names, normal.file.names)
      yaml.list <- setNames(as.list(all.files), all.file.names)
      # TODO : add the pairings
      pair.list <- setNames(as.list(normal.file.names), file.names)
      yaml.structure <- list(samples = yaml.list, pairings = pair.list)
    }
    write_yaml(yaml.structure, save.path)
}

compute.index <- \(path.to.bam = "ichorCNA-test/") {
  file.names <- list.files(path.to.bam, pattern = "\\.bam$")
  print(file.names)
  for (file in file.names) {
    bai.name <- paste0(file, ".bai")
    bai.exists <- file.exists(paste0(path.to.bam, bai.name))
    if (!bai.exists) {
      print(paste0(bai.name, " does not exists. Computing index.."))
      system(paste0("samtools index -b ", paste0(path.to.bam, file), 
                    paste0(" ", path.to.bam, file, ".bai")))
    }
  }
}

