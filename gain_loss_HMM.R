library(HMM)
library(data.table)
library(preprocessCore)


#assume a minimu size of 10, which is correct for the
#bigwigs we are using
#change this if the resolution is different
#by removing the aggregate we can speed up the analysis
read.chip.fast <- function( file, res = 10){
	#use data table to quickly read the bedgraph for an entire chromosome
	bdg <- fread(file, data.table=F)
	size <- bdg[,3]/res - bdg[,2]/res
	size[length(size)] <- size[length(size)] + 10-floor(sum(size))%%10
	score <- rep(bdg[,4], size)
	score.matrix <- matrix(score, byrow=T, ncol=10)
	cov.vec <- apply(score.matrix, 1, mean)
	#combine the coverage vector with a position vector, which are 100bp bins
	cov <- cbind(seq(0,length(cov.vec)-1, 1)*100, cov.vec)
	cov
}	

#use bigWigToBedGraph from the UCSC to generate a bedgraph for a specific chromosome
#note that the user needs to provide the path to the bigWigToBedGraph program on the system
#download this for your system here: http://hgdownload.soe.ucsc.edu/admin/exe/
read.bigwig <- function( bw.file, chrom ){
	command <- paste0("~/bin/bigWigToBedGraph ", bw.file, " temp.bdg -chrom=", chrom)
	system(command)
}	


#bigWigs from GSE135180
bw.0.file <- "GSM3992901_5_WaplC6-0h_antiRad21_sample_calibrated.bw"
bw.24.file <- "GSM3992902_6_WaplC6-24h_antiRad21_sample_calibrated.bw"


peak.bed <- NULL

#loop over chromosomes and determine gained, lost and non-changing regions for 
for( chr in paste0("chr", c(1:19,"X")) ){
	cat(chr, "\n")

	read.bigwig( bw.0.file, chr)
	ch0 <- read.chip.fast( "temp.bdg" )
	
	read.bigwig( bw.24.file, chr)
	ch24 <- read.chip.fast( "temp.bdg" )
	
	cat("Done reading\n")	
	#merge the two datasets and set missing values to 0
	ch <- merge(ch0, ch24, by=1, all=T)
	ch[!is.finite(ch[,2]),2] <- 0
	ch[!is.finite(ch[,3]),3] <- 0

	#perform quantile normalization
	ch[,2:3] <- normalize.quantiles(as.matrix(ch[,2:3]))

	States <- c("loss", "gain", "no_change")
	Symbols <- c("chip_up", "chip_down", "chip_same")
	
	#trans=0.001
	trans = 1/1e6

	
	#note that chip_up means a higher signal in 24h vs untreated (0h)

	#emission probs for the gain state
	gain.up.em <- 0.315116226
	gain.down.em <- 0.010701464
	gain.same.em <- 0.674182310

	#emission probs for the loss state
	loss.up.em <- 0.004151122
	loss.down.em <- 0.213600871
	loss.same.em <- 0.782248006

	#emission probs for the no change state
	nc.up.em <- 0.037
	nc.down.em <- 0.037
	nc.same.em <- 0.926

	#set a threshold for what is considered when the signal is up or down
	chip.th <- 1 
	
	#create the matrices with the transition probabilities and the emission probabilities
	transProbs <- matrix(c(1-trans, trans, trans, trans, 1-trans, trans, trans, trans, 1-trans ), c(length(States), length(States)), byrow = TRUE)
	emissionProbs <- matrix(c(loss.up.em, loss.down.em, loss.same.em,
														 gain.up.em, gain.down.em, gain.same.em,
														 nc.up.em, nc.down.em, nc.same.em), 
														 nrow=c(length(States), ncol=length(Symbols)), byrow = TRUE)
	
	
	#initialize the hidden markov model
	hmm <- initHMM( States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs)
	
	#calculate the difference of the 24h and untreated ChIP signal
	chip.diff <- ch[,3]-ch[,2]
	
	#initiate a vection in which all the chip signals are the same
	obs <- rep("chip_same", nrow(ch) )
	#set which bins are considered increasing and decreasing
	obs[chip.diff > chip.th] <- "chip_up"
	obs[chip.diff < -chip.th] <- "chip_down"
	
	#perform the viterbi decoding based on the observations
	vit <- viterbi(hmm, obs)
	
	#get the gain domains
	xv <- which( vit == 'gain' )
	pos <- ch[xv,1]
	i <- which(diff(pos) > 100)
	
	#get the positions in the genome
	start <- c(pos[1],pos[i+1])
	end <- c(pos[i],pos[length(pos)])

	#create bed formatted data frame
	chrom.bed <- data.frame(chrom=chr, start=start, end=end, type="gain")
	
	peak.bed <- rbind(peak.bed, chrom.bed)
	
	#get the loss domains
	xv <- which( vit == 'loss' )
	pos <- ch[xv,1]
	i <- which(diff(pos) > 100)
	
	#get the positions in the genome
	start <- c(pos[1],pos[i+1])
	end <- c(pos[i],pos[length(pos)])
	
	#create bed formatted data frame
	chrom.bed <- data.frame(chrom=chr, start=start, end=end, type="loss")
	
	#concatenate gained and lost domains into on data.frame
	peak.bed <- rbind(peak.bed, chrom.bed)

}

#write results to bed file
#write.table(peak.bed, "gain_loss_out.bed", col=F, row=F, quote=F, sep="\t")
