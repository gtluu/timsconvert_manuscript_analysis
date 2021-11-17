library(Cardinal)
library(devtools)
library(cardinalscripts)
library(stringr)
library(fuzzyjoin)
library(dplyr)
library(ggplot2)
library(reticulate)
kneed <- import('kneed')
numpy <- import('numpy')

# Set number of workers to 2.
register(SnowParam(workers=2, log=TRUE, logdir=getwd()))

# Load in imzML format data from "data" folder.
rugose_tsf <- readMSIData('20210920_vc_rugose_1mMTCA_gordon.imzML',
                          resolution=200,
                          units='ppm')
rugose_tsf_spots <- '20210920_vc_rugose_1mMTCA_gordon_spot_list.txt'
rugose_tsf_rois <- '20210920_vc_rugose_1mMTCA_gordon.csv'

# Update metadata in pixelData() via cardinalscripts::updateMetadata().
rugose_tsf <- updateMetadata(rugose_tsf,
                             spotFile=rugose_tsf_spots,
                             roiFile=rugose_tsf_rois)
centroided(rugose_tsf) <- FALSE

# Preprocessing
# Normalization, Signal Smoothing, and Baseline Subtraction
preprocessed <- normalize(rugose_tsf, method='tic')
preprocessed <- smoothSignal(preprocessed, method='sgolay')
preprocessed <- reduceBaseline(preprocessed, method='median')
preprocessed <- process(preprocessed)
# Peak Picking, Alignment, and Filtering
peaks <- peakPick(preprocessed, method='simple', SNR=3)
peaks <- peakAlign(peaks, tolerance=NA, units='mz')
peaks <- peakFilter(peaks, freq.min=0.05)
processed_rugose_tsf <- process(peaks)

# Generate multivariate segmentation models via spatial shrunken centroids with
# varying r, k, and s parameters.
set.seed(1995)

rparam <- c(1,2,3)
kparam <- c(2,4,6,8,10,12,14,16,18,20)
sparam <- c(0,3,6,9,12,15)

rugose_tsf_ssc <- spatialShrunkenCentroids(processed_rugose_tsf,
                                           r=rparam,
                                           k=kparam,
                                           s=sparam)

# Optimal r, k, and s parameters were chosen after manual inspection of SSC
# models.
rugose_tsf_ssc_params$params$r <- 2
rugose_tsf_ssc_params$params$k <- 4
rugose_tsf_ssc_params$params$s <- 3

# Get table of t-statistics for SSC model where r=2, k=4, and s=3 to identify
# signals contributing to each segment.
rugose_tsf_ssc_model <- modelData(rugose_tsf_ssc)
rownames(rugose_tsf_ssc_model) <- c(1:nrow(rugose_tsf_ssc_model))
index <- as.numeric(rownames(rugose_tsf_ssc_model)[rugose_tsf_ssc_model$r == 2 &
                                                   rugose_tsf_ssc_model$k == 4 &
                                                   rugose_tsf_ssc_model$s == 3])
rugose_tsf_ssc_stat_table <- as.data.frame(resultData(rugose_tsf_ssc)[[index]]$statistic)

# Export SSC images for all models to PDF.
pdf('20210920_vc_rugose_1mMTCA_gordon.pdf')
for (r in rparam) {
  for (k in kparam) {
    for (s in sparam) {
      print(image(rugose_tsf_ssc, model=list(r=r, k=k, s=s)))
    }
  }
}
dev.off()

# Output images used in Figure S3C.
image(rugose_tsf_ssc, model=list(r=2, k=4, s=3), xlim=c(10, 45))
