# normalization bigWig function
# BWdir0h (directory of input wildtype bigWig)
# BWdir24h (directory of input mutant bigWig)
# BW0h (file name of input wildtype bigWig)
# BW24h (file name of input mutant bigWig)
# refdir (directory of deeptools alignment for all the peaks from internal control)
# refmatrix (matrix of deeptools alignment for all the peaks from internal control)
# outdir (output directory)
# out0h (output file name of normalized wildtype bigWig)
# out24h (output file name of normalized mutant bigWig)

SICHIP_norm <- function(BWdir0h, BWdir24h, BW0h, BW24h, refdir, refmatrix, outdir, out0h, out24h){
  # load the required libraries
  library(rtracklayer)
  
  # load bigwig files
  bigwig0h <- import.bw(paste0(BWdir0h, BW0h))
  bigwig24h <- import.bw(paste0(BWdir24h, BW24h))
  
  # load reference matrix calculated with the same method (reference is the Hek293 spike-in)
  ref <- data.frame(read.table(paste0(refdir, refmatrix), skip=1, header=F, sep="\t"))
  ref_ratio <- rbind(apply(ref[, 7:306], 2, mean)/apply(ref[, 7:306], 2, mean), 
                     apply(ref[, 7:306], 2, mean)/apply(ref[, 307:606], 2, mean)
  )
  row.names(ref_ratio) <- c("w0h", "w24h")
  ref_ratio <- cbind(apply(ref_ratio, 1, function(x)mean(x[which(!is.na(x))])), ref_ratio)
  
  # normalize data matrix by spike-in reference
  bigwig0h$score <- bigwig0h$score*ref_ratio[1,1]
  bigwig24h$score <- bigwig24h$score*ref_ratio[2,1]
  
  # export normalized bigWig files
  export.bw(bigwig0h, paste0(outdir,out0h))
  export.bw(bigwig24h, paste0(outdir,out24h))
}

## Example (mm10):
# SICHIP_norm("/DATA/projects/ES_Wapl_AID/ChIPseq/WCR_ChIP/bigWig_sample/",
#             "/DATA/projects/ES_Wapl_AID/ChIPseq/WCR_ChIP/bigWig_sample/",
#             "1_Wapl-0D-antiWapl_SF.3741_MQ15_sample.bw",
#             "4_Wapl-1D-antiWapl_SF.2194_MQ15_sample.bw",
#             "/DATA/projects/ES_Wapl_AID/ChIPseq/WCR_ChIP/Peak_reference/",
#             "Hek_Wapl-0D_1D-antiWapl_reference_scale.mat.gz",
#             "/DATA/projects/ES_Wapl_AID/ChIPseq/Normalized_SIChIP/",
#             "1_WaplC6-0h_antiWapl_sample.bw",
#             "1_WaplC6-24h_antiWapl_sample.bw")
