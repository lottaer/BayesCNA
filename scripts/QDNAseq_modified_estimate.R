library(methods)
library(QDNAseq)
library(ggplot2)
library(magrittr)
library(rstudioapi)

#Set working directory
getActiveDocumentContext()$path %>% dirname %>% setwd %T>% print 

#Set subfolders for study to loop through
sampleDir <- 'cell_line_subs'
cat('Performing analysis on samples in folder ',sampleDir)

#Define directories: for intermediate plotting, input and outputs
plotdir <- sampleDir
inputdir <- sampleDir
outputdir <- sampleDir

#Read in all bams in input directory in one go, they will be stored as separate samples in the analysis
#and all steps will be executed to each one of them
bins <- getBinAnnotations(binSize = 500)
readCounts <- binReadCounts(bins, path = inputdir)
# Manually define chromosomes that will be used for estimating correction
chromList <- c('4','9','15','16')

#Filter profiles and plot isobar plots to see dependence of median read  count on GC content and mappability	
readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE, mappability=75)

# Estimate GC content correction, and plot relationship between observed standard deviation and read depth
# low-quality DNA will be noisier than the linear line and will appear above it
# The object created here is then modified if estimation is carried out "manually"
readCountsFiltered <- estimateCorrection(readCountsFiltered)

# Do modified estimate fitting based on selected (mostly normal) chromosomes alone
x <- as.data.frame(readCountsFiltered@assayData$counts)
samples <- names(x)
x$chrom <- sapply(rownames(x), function(z) unlist(strsplit(z,':'))[1])
x[,c('gc','mappability','use')] <- readCountsFiltered@featureData@data[match(rownames(x), rownames(readCountsFiltered@featureData@data)),c('gc','mappability','use')]
x$use <- x$use & !is.na(x$mappability) & !is.na(x$gc) & (x$gc>25) & (x$gc<90)
x$gc <- round(x$gc)
x$mappability <- round(x$mappability)
all.combinations <- expand.grid(gc=unique(x$gc[!is.na(x$gc)]),
                                mapp=unique(x$mappability[!is.na(x$mappability)]))
rownames(all.combinations) <- paste0(all.combinations$gc, "-",
                                     all.combinations$mapp)

y <- subset(x, chrom %in% chromList)

samples <- samples[-(which(samples == "E5_2-10"))]

for(var in samples){
  xx <- data.frame(gc=y[y$use,'gc'], map=y[y$use,'mappability'], counts=y[y$use,var])
  median.counts <- aggregate(xx[,'counts'], by=list(gc=xx$gc, mapp=xx$map), FUN=median)
  l <- loess(x~gc*mapp, data=median.counts, span=0.8, family="symmetric")

  # Fill in fit for all elements
  fit <- as.vector(predict(l, all.combinations)); names(fit) <- rownames(all.combinations)
  x[,paste0(var,'_fit')] <- fit[paste0(x$gc,'-',x$mappability)]

  # Overwrite the corresponding object within the QDNAseq pipeline
  unlockBinding('fit',readCountsFiltered@assayData)
  readCountsFiltered@assayData$fit[,var] <- x[,paste0(var,'_fit')]
  lockBinding('fit',readCountsFiltered@assayData)
}

save(readCountsFiltered,file = paste0(outputdir,"/Read_counts_raw_500.RData",sep=''))
#Correct bins using the estimator and normalise and smooth the aquired copy numbers
#Plot copy number profiles
copyNumbers <- correctBins(readCountsFiltered)

copyNumbersNormalized <- normalizeBins(copyNumbers, method='mode') # mode-based normalisation works better for this data
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
pdf(paste0(plotdir,"/QDNAseq_copynumbers_smooth.pdf",sep=''),5.75,4)
plot(copyNumbersSmooth)
dev.off()

#Perform segmentation and normalise segment, then plot the segmented copy numbers
copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
#copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
pdf(paste0(plotdir,"/QDNAseq_copynumbers_segments.pdf",sep=''),5.75,4)
plot(copyNumbersSegmented)
dev.off()
save(copyNumbersSegmented,file=paste0(outputdir,"/Copynumber_calls_500.RData"))


