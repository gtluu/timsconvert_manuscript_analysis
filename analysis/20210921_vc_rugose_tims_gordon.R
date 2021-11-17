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
rugose_tdf <- readMSIData('data/20210921_vc_rugose_tims_gordon.imzML',
                          resolution=200,
                          units='ppm')
rugose_tdf_spots <- 'data/20210921_vc_rugose_tims_gordon_spot_list.txt'
rugose_tdf_rois <- 'data/20210921_vc_rugose_tims_gordon.csv'

# Update metadata in pixelData() via cardinalscripts::updateMetadata().
rugose_tdf <- updateMetadata(rugose_tdf,
                             spotFile=rugose_tdf_spots,
                             roiFile=rugose_tdf_rois)
centroided(rugose_tdf) <- FALSE

# Preprocessing
# Normalization, Signal Smoothing, and Baseline Subtraction
preprocessed <- normalize(rugose_tdf, method='tic')
preprocessed <- smoothSignal(preprocessed, method='sgolay')
preprocessed <- reduceBaseline(preprocessed, method='median')
preprocessed <- process(preprocessed)
# Peak Picking, Alignment, and Filtering
peaks <- peakPick(preprocessed, method='mad', SNR=3)
peaks <- peakAlign(peaks, tolerance=NA, units='mz')
peaks <- peakFilter(peaks, freq.min=0.05)
processed_rugose_tdf <- process(peaks)

# Generate multivariate segmentation models via spatial shrunken centroids with
# varying r, k, and s parameters.
set.seed(1995)

rparam <- c(1,2,3)
kparam <- c(2,4,6,8,10,12,14,16,18,20)
sparam <- c(0,3,6,9,12,15)

rugose_tdf_ssc <- spatialShrunkenCentroids(processed_rugose_tdf,
                                           r=rparam,
                                           k=kparam,
                                           s=sparam)

# Optimal r, k, and s parameters were chosen after manual inspection of SSC
# models.
rugose_tdf_ssc_params$params$r <- 2
rugose_tdf_ssc_params$params$k <- 6
rugose_tdf_ssc_params$params$s <- 0

# Get table of t-statistics for SSC model where r=2, k=6, and s=0 to identify
# signals contributing to each segment.
rugose_tdf_ssc_model <- modelData(rugose_tdf_ssc)
rownames(rugose_tdf_ssc_model) <- c(1:nrow(rugose_tdf_ssc_model))
index <- as.numeric(rownames(rugose_tdf_ssc_model)[rugose_tdf_ssc_model$r == 2 &
                                                   rugose_tdf_ssc_model$k == 6 &
                                                   rugose_tdf_ssc_model$s == 0])
rugose_tdf_ssc_stat_table <- as.data.frame(resultData(rugose_tdf_ssc)[[index]]$statistic)

# Export SSC images for all models to PDF.
pdf('20210921_vc_rugose_tims_gordon.pdf')
for (r in rparam) {
  for (k in kparam) {
    for (s in sparam) {
      print(image(rugose_tdf_ssc, model=list(r=r, k=k, s=s)))
    }
  }
}
dev.off()

# Output images used in Figure S3A and S3B.
image(rugose_tdf_ssc, model=list(r=2, k=6, s=0), xlim=c(7, 40))
image(processed_rugose_tdf, mz=mz(processed_rugose_tdf)[462], xlim=c(7, 40))
