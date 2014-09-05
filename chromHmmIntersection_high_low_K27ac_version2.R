cols <- c("purple4","red","blue","green","yellow","purple","orange",
          "pink","black","grey","rosybrown","lightblue",
          "lightgreen","hotpink","brown","wheat3","snow",
          "maroon1","darkblue","darkgreen","tomato",
          "turquoise","orchid","orangered4","yellowgreen","magenta")
cell <- commandArgs(TRUE)[1]
pdf(paste("chromHmm.",cell,".V2.pdf",sep=""))
par(xpd=T, mar=par()$mar+c(0,0,0,10))

random_distal_H3K27ac_table <- as.data.frame(read.table(paste("mean_H3K27ac.",cell,".randombg3.distal.txt",sep="")))
dnase_distal_H3K27ac_table <- as.data.frame(read.table(paste("mean_H3K27ac.",cell,".distal.txt",sep="")))
total_peaks_distal_dnase_peaks <- nrow(dnase_distal_H3K27ac_table)
total_peaks_in_distal_random_peaks <- nrow(random_distal_H3K27ac_table)
cat(total_peaks_distal_dnase_peaks == total_peaks_in_distal_random_peaks,"\n")

random_proximal_H3K27ac_table <- as.data.frame(read.table(paste("mean_H3K27ac.",cell,".randombg3.proximal.txt",sep="")))
dnase_proximal_H3K27ac_table <- as.data.frame(read.table(paste("mean_H3K27ac.",cell,".proximal.txt",sep="")))
total_peaks_in_proximal_dnase_peaks <- nrow(dnase_proximal_H3K27ac_table)
total_peaks_in_proximal_random_peaks <- nrow(random_proximal_H3K27ac_table)
cat(total_peaks_in_proximal_dnase_peaks == total_peaks_in_proximal_random_peaks,"\n")

bg_H3K27_distal_values <- as.numeric(random_distal_H3K27ac_table[,ncol(random_distal_H3K27ac_table)])
iles_5 <- quantile(bg_H3K27_distal_values, prob=seq(0,1,0.05))
distal_threshold <- as.numeric(iles_5[length(iles_5)-1]) # last ile

bg_H3K27_proximal_values <- as.numeric(random_proximal_H3K27ac_table[,ncol(random_proximal_H3K27ac_table)])
iles_5 <- quantile(bg_H3K27_proximal_values, prob=seq(0,1,0.05))
proximal_threshold <- as.numeric(iles_5[length(iles_5)-1]) # last ile
		
distal_dnase_high <- as.data.frame(read.table(paste(cell,"_dnase_high_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_dnase_high[,2] <- distal_dnase_high[,2]/total_peaks_distal_dnase_peaks
distal_dnase_high <- distal_dnase_high[order(distal_dnase_high[,1]),]

distal_dnase_low <- as.data.frame(read.table(paste(cell,"_dnase_low_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_dnase_low[,2] <- distal_dnase_low[,2]/total_peaks_distal_dnase_peaks
distal_dnase_low <- distal_dnase_low[order(distal_dnase_high[,1]),]

distal_random_high <- as.data.frame(read.table(paste(cell,"_random_high_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_random_high[,2] <- distal_random_high[,2]/total_peaks_in_distal_random_peaks
distal_random_high <- distal_random_high[order(distal_dnase_high[,1]),]
distal_random_low <- as.data.frame(read.table(paste(cell,"_random_low_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_random_low[,2] <- distal_random_low[,2]/total_peaks_in_distal_random_peaks
distal_random_low <- distal_random_low[order(distal_dnase_high[,1]),]

distal <- cbind(distal_dnase_high[,2],distal_random_high[,2],distal_dnase_low[,2],distal_random_low[,2])
rownames(distal) <- distal_dnase_high[,1]
bp_coords <- barplot(distal, beside=FALSE, cex.main=0.9, col=cols, las=2, main=paste(cell, "distal"), xaxt="n", ylim=c(0,1), cex.axis=0.6)
axis(1, at=bp_coords, c("dnase(top 5% H3K27ac)", "dnase(other)", "rand(top 5% H3K27ac)","rand(other)"), cex=0.6, cex.axis=0.6)
legend(5,1, rownames(distal), col=cols,pch=15, bty="n", cex=0.6)



proximal_dnase_high <- as.data.frame(read.table(paste(cell,"_dnase_high_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))	
proximal_dnase_high[,2] <- proximal_dnase_high[,2]/total_peaks_in_proximal_dnase_peaks
proximal_dnase_high <- proximal_dnase_high[order(proximal_dnase_high[,1]),]

proximal_dnase_low <- as.data.frame(read.table(paste(cell,"_dnase_low_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))	
proximal_dnase_low[,2] <- proximal_dnase_low[,2]/total_peaks_in_proximal_dnase_peaks
proximal_dnase_low <- proximal_dnase_low[order(proximal_dnase_high[,1]),]

proximal_random_high <- as.data.frame(read.table(paste(cell,"_random_high_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))	
proximal_random_high[,2] <- proximal_random_high[,2]/total_peaks_in_proximal_random_peaks
proximal_random_high <- proximal_random_high[order(proximal_dnase_high[,1]),]

proximal_random_low <- as.data.frame(read.table(paste(cell,"_random_low_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))	
proximal_random_low[,2] <- proximal_random_low[,2]/total_peaks_in_proximal_random_peaks
proximal_random_low <- proximal_random_low[order(proximal_dnase_high[,1]),]

proximal <- cbind(proximal_dnase_high[,2],proximal_random_high[,2],proximal_dnase_low[,2],proximal_random_low[,2])
rownames(proximal) <- proximal_dnase_high[,1]
bp_coords <- barplot(proximal, beside=FALSE, cex.main=0.9, col=cols, las=2, main=paste(cell, "proximal"), xaxt="n", ylim=c(0,1), cex.axis=0.6)
axis(1, at=bp_coords, c("dnase(top 5% H3K27ac)", "dnase(other)", "rand(top 5% H3K27ac)","rand(other)"), cex=0.6, cex.axis=0.6)
legend(5,1, rownames(proximal), col=cols,pch=15, bty="n", cex=0.6)

dev.off()
