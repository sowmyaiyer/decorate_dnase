cells <- scan("cells_with_H3K27Ac_and_Dnase.test.txt", what="character")
cell <- commandArgs(TRUE)[1]
pdf(paste("mean_H3K27Ac_matched_random_",cell,".pdf",sep=""))
#pdf(paste("mean_H3K27Ac_matched_random.pdf",sep=""))
#for (cell in cells)
{
	cat(cell,"\n")
	par(mfrow=c(2,2))
	dist_random <- as.matrix(read.table(paste("mean_H3K27ac.",cell,".randombg3.distal.txt",sep="")))
	prox_random <- as.matrix(read.table(paste("mean_H3K27ac.",cell,".randombg3.proximal.txt",sep="")))
	
	dist_dnase <- as.matrix(read.table(paste("mean_H3K27ac.",cell,".distal.txt",sep="")))
	prox_dnase <- as.matrix(read.table(paste("mean_H3K27ac.",cell,".proximal.txt",sep="")))
	ymax <- length(c(as.numeric(dist_dnase[,5]), as.numeric(dist_random[,5])))
	xmax <- max(c(as.numeric(dist_dnase[,5])), as.numeric(dist_random[,5]))	
	dnase_hist_dist <- hist(as.numeric(dist_dnase[,5]), breaks=seq(0,500000,1), plot=FALSE)
	random_hist_dist <- hist(as.numeric(dist_random[,5]),  breaks=seq(0,500000,1), plot=FALSE)


	sum_last_bin_dist_dnase <- sum(dnase_hist_dist$counts[12:length(dnase_hist_dist$counts)])
	sum_last_bin_dist_random <- sum(random_hist_dist$counts[12:length(random_hist_dist$counts)])
	break_last_bin_dist <- ">= 11"
	counts_dist_dnase <- c(dnase_hist_dist$counts[1:11],sum_last_bin_dist_dnase)
	counts_dist_random <- c(random_hist_dist$counts[1:11],sum_last_bin_dist_random)
	breaks_dist <- c(dnase_hist_dist$breaks[1:11],break_last_bin_dist)

	bp_coords <- barplot(rbind(counts_dist_dnase,counts_dist_random), col=c("darkgreen","lightgreen"), beside=TRUE, axes=FALSE, main=paste(cell,"distal"), cex.main=0.7,space=c(0.3,1), border=NA, xlab="H3K27Ac", cex.lab=0.75)
	axis(1, at=colMeans(bp_coords), breaks_dist, cex.axis=0.60, las=1, tck=FALSE, lwd=0, line=-1, las=2)
        axis(2, cex=0.60, cex.axis=0.75, las=1, col="black", col.axis="black")
	legend("topright", c(paste("random n =",nrow(dist_random)), paste("dnase n =",nrow(dist_dnase))), pch=15, col=c("lightgreen","darkgreen","white","white"), bty="n", cex=0.75)

	dnase_hist_prox <- hist(as.numeric(prox_dnase[,5]), breaks=seq(0,500000,1), plot=FALSE )
        random_hist_prox <- hist(as.numeric(prox_random[,5]),  breaks=seq(0,500000,1), plot=FALSE)

        sum_last_bin_prox_dnase <- sum(dnase_hist_prox$counts[12:length(dnase_hist_prox$counts)])
        sum_last_bin_prox_random <- sum(random_hist_prox$counts[12:length(random_hist_prox$counts)])
        break_last_bin_prox <- " >= 21"
        counts_prox_dnase <- c(dnase_hist_prox$counts[1:21],sum_last_bin_prox_dnase)
        counts_prox_random <- c(random_hist_prox$counts[1:21],sum_last_bin_prox_random)
	breaks_prox <- c(dnase_hist_prox$breaks[1:21],break_last_bin_prox)

        bp_coords <- barplot(rbind(counts_prox_dnase,counts_prox_random), col=c("darkblue","lightblue"), beside=TRUE, axes=FALSE, main=paste(cell,"proximal"), cex.main=0.7,space=c(0.3,1), border=NA, xlab="H3K27Ac", cex.lab=0.75)
        axis(1, at=colMeans(bp_coords), breaks_prox, cex.axis=0.60, las=1, tck=FALSE, lwd=0, line=-1, las=2)
        axis(2, cex=0.60, cex.axis=0.75, las=1, col="black", col.axis="black")
	legend("topright", c(paste("random n =",nrow(prox_random)), paste("dnase n =",nrow(prox_dnase))), pch=15, col=c("lightblue","darkblue","white","white"), bty="n", cex=0.75)

}
dev.off()
