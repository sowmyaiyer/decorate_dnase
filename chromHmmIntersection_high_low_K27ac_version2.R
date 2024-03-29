#cols <- c("gray","gray","rosybrown","rosybrown","yellow","yellow","pink",
#         "pink","lightblue","lightblue","lightblue","lightblue",
#          "yellow","pink","yellowgreen","pink","gray",
#          "yellowgreen","yellowgreen","yellowgreen","darkslategray4",
#          "darkslategray4","darkslategray4","darkslategray4","yellowgreen","yellowgreen")

grouping <- list("other"=c(".","Art","Low"), 
		 "Ctcf"=c("Ctcf", "CtcfO"), 
		 "Dnase/Faire"=c("DnaseD", "DnaseU","FaireW"), 
		 "Elon/H4K20"=c("Elon", "ElonW","H4K20","Gen3'"),
		 "enhancer"=c("Enh","EnhF","EnhW","EnhWF"),
		 "promoter"=c("TssF","Tss","PromP","PromF","Pol2","Gen5'"),
		 "quiescent"=c("Quies"),
		 "repressed"=c("Repr","ReprD","ReprW"))

cols_list <- list("other"="gray", 
		  "Ctcf"="rosybrown", 
		  "Dnase/Faire"="yellow", 
		  "Elon/H4K20"="pink",
		  "enhancer"="lightblue",
		  "promoter"="yellowgreen",
		  "quiescent"="lightslategray",
		  "repressed"="darkslategray4")


cell <- commandArgs(TRUE)[1]
pdf(paste("chromHmm.",cell,".withChipSeq.pdf",sep=""))
par(xpd=T, mar=par()$mar+c(0,0,0,15))

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


mat_barplot_distal <- matrix(nrow=length(grouping), ncol=4)
rownames(mat_barplot_distal) <- names(grouping)
colnames(mat_barplot_distal) <- c("dnase_high_H3K27ac","dnase_low_H3K27ac", "random_high_H3K27ac","random_low_H3K27ac")

mat_barplot_chipseq_distal <- matrix(nrow=length(grouping), ncol=4)
rownames(mat_barplot_chipseq_distal) <- names(grouping)
colnames(mat_barplot_chipseq_distal) <- c("dnase_high_H3K27ac_chipseq","dnase_low_H3K27ac_chipseq", "random_high_H3K27ac_chipseq","random_low_H3K27ac_chipseq")

distal_dnase_high <- as.data.frame(read.table(paste(cell,"_dnase_high_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_dnase_high[,2] <- distal_dnase_high[,2]/total_peaks_distal_dnase_peaks
distal_dnase_high <- distal_dnase_high[order(distal_dnase_high[,1]),]

distal_dnase_high_chipseq <- as.data.frame(read.table(paste(cell,"_dnase_high_H3K27ac_intersect_chipSeq_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_dnase_high_chipseq[,2] <- distal_dnase_high_chipseq[,2]/total_peaks_distal_dnase_peaks
distal_dnase_high_chipseq <- distal_dnase_high_chipseq[order(distal_dnase_high_chipseq[,1]),]
for (group in rownames(mat_barplot_distal))
{
	mat_barplot_distal[group,"dnase_high_H3K27ac"] <- sum(distal_dnase_high[match(grouping[[group]],distal_dnase_high[,1]),2])
	mat_barplot_chipseq_distal[group,"dnase_high_H3K27ac_chipseq"] <- sum(distal_dnase_high_chipseq[match(grouping[[group]],distal_dnase_high_chipseq[,1]),2])
}

distal_dnase_low <- as.data.frame(read.table(paste(cell,"_dnase_low_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_dnase_low[,2] <- distal_dnase_low[,2]/total_peaks_distal_dnase_peaks
distal_dnase_low <- distal_dnase_low[order(distal_dnase_high[,1]),]

distal_dnase_low_chipseq <- as.data.frame(read.table(paste(cell,"_dnase_low_H3K27ac_intersect_chipSeq_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_dnase_low_chipseq[,2] <- distal_dnase_low_chipseq[,2]/total_peaks_distal_dnase_peaks
distal_dnase_low_chipseq <- distal_dnase_low_chipseq[order(distal_dnase_high_chipseq[,1]),]
for (group in rownames(mat_barplot_distal))
{
        mat_barplot_distal[group,"dnase_low_H3K27ac"] <- sum(distal_dnase_low[match(grouping[[group]],distal_dnase_low[,1]),2])
        mat_barplot_chipseq_distal[group,"dnase_low_H3K27ac_chipseq"] <- sum(distal_dnase_low_chipseq[match(grouping[[group]],distal_dnase_low_chipseq[,1]),2])
}

distal_random_high <- as.data.frame(read.table(paste(cell,"_random_high_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_random_high[,2] <- distal_random_high[,2]/total_peaks_in_distal_random_peaks
distal_random_high <- distal_random_high[order(distal_dnase_high[,1]),]

distal_random_high_chipseq <- as.data.frame(read.table(paste(cell,"_random_high_H3K27ac_intersect_chipSeq_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_random_high_chipseq[,2] <- distal_random_high_chipseq[,2]/total_peaks_in_distal_random_peaks
distal_random_high_chipseq <- distal_random_high_chipseq[order(distal_dnase_high_chipseq[,1]),]
for (group in rownames(mat_barplot_distal))
{
        mat_barplot_distal[group,"random_high_H3K27ac"] <- sum(distal_random_high[match(grouping[[group]],distal_random_high[,1]),2])
        mat_barplot_chipseq_distal[group,"random_high_H3K27ac_chipseq"] <- sum(distal_random_high_chipseq[match(grouping[[group]],distal_random_high_chipseq[,1]),2])
}

distal_random_low <- as.data.frame(read.table(paste(cell,"_random_low_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_random_low[,2] <- distal_random_low[,2]/total_peaks_in_distal_random_peaks
distal_random_low <- distal_random_low[order(distal_dnase_high[,1]),]

distal_random_low_chipseq <- as.data.frame(read.table(paste(cell,"_random_low_H3K27ac_intersect_chipSeq_intersect_chromHmm.25.states.distal.counts.txt",sep="")))
distal_random_low_chipseq[,2] <- distal_random_low_chipseq[,2]/total_peaks_in_distal_random_peaks
distal_random_low_chipseq <- distal_random_low_chipseq[order(distal_dnase_high_chipseq[,1]),]
for (group in rownames(mat_barplot_distal))
{
        mat_barplot_distal[group,"random_low_H3K27ac"] <- sum(distal_random_low[match(grouping[[group]],distal_random_low[,1]),2])
        mat_barplot_chipseq_distal[group,"random_low_H3K27ac_chipseq"] <- sum(distal_random_low_chipseq[match(grouping[[group]],distal_random_low_chipseq[,1]),2])
}
bp_coords <- barplot(mat_barplot_distal, beside=FALSE, cex.main=0.9, col=unlist(cols_list), las=2, main=paste(cell, "distal"), xaxt="n", ylim=c(0,1), cex.axis=0.6, border=NA, space=c(2,2,2,2))
barplot(as.matrix(mat_barplot_chipseq_distal[,"dnase_high_H3K27ac_chipseq"]),  space=bp_coords[1]+0.5, beside=FALSE, cex.main=0.9, col=unlist(cols_list), axes=FALSE, add=TRUE, border=NA)
barplot(as.matrix(mat_barplot_chipseq_distal[,"dnase_low_H3K27ac_chipseq"]),  space=bp_coords[2]+0.5, beside=FALSE, cex.main=0.9, col=unlist(cols_list), axes=FALSE, add=TRUE, border=NA)
barplot(as.matrix(mat_barplot_chipseq_distal[,"random_high_H3K27ac_chipseq"]),  space=bp_coords[3]+0.5, beside=FALSE, cex.main=0.9, col=unlist(cols_list), axes=FALSE, add=TRUE, border=NA)
barplot(as.matrix(mat_barplot_chipseq_distal[,"random_low_H3K27ac_chipseq"]),  space=bp_coords[4]+0.5, beside=FALSE, cex.main=0.9, col=unlist(cols_list), axes=FALSE, add=TRUE, border=NA)
#legend(15,0.8, c(names(cols_list),"fraction with ChIP-seq intersection"), col=c(unlist(cols_list),"black"), pch=15, bty="n", cex=0.8, pt.cex=1.2)
text(x=bp_coords+0.5,  par("usr")[3]-0.01, labels=colnames(mat_barplot_distal), cex=0.6, srt = 30, pos = 2, xpd = TRUE)
legend(15,0.8, names(cols_list), col=unlist(cols_list), pch=15, bty="n", cex=0.8, pt.cex=1.2)


mat_barplot_proximal <- matrix(nrow=length(grouping), ncol=4)
rownames(mat_barplot_proximal) <- names(grouping)
colnames(mat_barplot_proximal) <-  c("dnase_high_H3K27ac","dnase_low_H3K27ac", "random_high_H3K27ac","random_low_H3K27ac")

mat_barplot_chipseq_proximal <- matrix(nrow=length(grouping), ncol=4)
rownames(mat_barplot_chipseq_proximal) <- names(grouping)
colnames(mat_barplot_chipseq_proximal) <- c("dnase_high_H3K27ac_chipseq","dnase_low_H3K27ac_chipseq", "random_high_H3K27ac_chipseq","random_low_H3K27ac_chipseq")


proximal_dnase_high <- as.data.frame(read.table(paste(cell,"_dnase_high_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))	
proximal_dnase_high[,2] <- proximal_dnase_high[,2]/total_peaks_in_proximal_dnase_peaks
proximal_dnase_high <- proximal_dnase_high[order(proximal_dnase_high[,1]),]

proximal_dnase_high_chipseq <- as.data.frame(read.table(paste(cell,"_dnase_high_H3K27ac_intersect_chipSeq_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))
proximal_dnase_high_chipseq[,2] <- proximal_dnase_high_chipseq[,2]/total_peaks_in_proximal_dnase_peaks
proximal_dnase_high_chipseq <- proximal_dnase_high_chipseq[order(proximal_dnase_high_chipseq[,1]),]
for (group in rownames(mat_barplot_proximal))
{
        mat_barplot_proximal[group,"dnase_high_H3K27ac"] <- sum(proximal_dnase_high[match(grouping[[group]],proximal_dnase_high[,1]),2])
        mat_barplot_chipseq_proximal[group,"dnase_high_H3K27ac_chipseq"] <- sum(proximal_dnase_high_chipseq[match(grouping[[group]],proximal_dnase_high_chipseq[,1]),2])
}

proximal_dnase_low <- as.data.frame(read.table(paste(cell,"_dnase_low_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))	
proximal_dnase_low[,2] <- proximal_dnase_low[,2]/total_peaks_in_proximal_dnase_peaks
proximal_dnase_low <- proximal_dnase_low[order(proximal_dnase_high[,1]),]

proximal_dnase_low_chipseq <- as.data.frame(read.table(paste(cell,"_dnase_low_H3K27ac_intersect_chipSeq_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))
proximal_dnase_low_chipseq[,2] <- proximal_dnase_low_chipseq[,2]/total_peaks_in_proximal_dnase_peaks
proximal_dnase_low_chipseq <- proximal_dnase_low_chipseq[order(proximal_dnase_high_chipseq[,1]),]
for (group in rownames(mat_barplot_proximal))
{
        mat_barplot_proximal[group,"dnase_low_H3K27ac"] <- sum(proximal_dnase_low[match(grouping[[group]],proximal_dnase_low[,1]),2])
        mat_barplot_chipseq_proximal[group,"dnase_low_H3K27ac_chipseq"] <- sum(proximal_dnase_low_chipseq[match(grouping[[group]],proximal_dnase_low[,1]),2])
}

proximal_random_high <- as.data.frame(read.table(paste(cell,"_random_high_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))	
proximal_random_high[,2] <- proximal_random_high[,2]/total_peaks_in_proximal_random_peaks
proximal_random_high <- proximal_random_high[order(proximal_dnase_high[,1]),]

proximal_random_high_chipseq <- as.data.frame(read.table(paste(cell,"_random_high_H3K27ac_intersect_chipSeq_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))
proximal_random_high_chipseq[,2] <- proximal_random_high_chipseq[,2]/total_peaks_in_proximal_random_peaks
proximal_random_high_chipseq <- proximal_random_high_chipseq[order(proximal_dnase_high_chipseq[,1]),]
for (group in rownames(mat_barplot_proximal))
{
        mat_barplot_proximal[group,"random_high_H3K27ac"] <- sum(proximal_random_high[match(grouping[[group]],proximal_random_high[,1]),2])
        mat_barplot_chipseq_proximal[group,"random_high_H3K27ac_chipseq"] <- sum(proximal_random_high_chipseq[match(grouping[[group]],proximal_random_high_chipseq[,1]),2])
}

proximal_random_low <- as.data.frame(read.table(paste(cell,"_random_low_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))	
proximal_random_low[,2] <- proximal_random_low[,2]/total_peaks_in_proximal_random_peaks
proximal_random_low <- proximal_random_low[order(proximal_dnase_high[,1]),]

proximal_random_low_chipseq <- as.data.frame(read.table(paste(cell,"_random_low_H3K27ac_intersect_chipSeq_intersect_chromHmm.25.states.proximal.counts.txt",sep="")))
proximal_random_low_chipseq[,2] <- proximal_random_low_chipseq[,2]/total_peaks_in_proximal_random_peaks
proximal_random_low_chipseq <- proximal_random_low_chipseq[order(proximal_dnase_high_chipseq[,1]),]
for (group in rownames(mat_barplot_proximal))
{
        mat_barplot_proximal[group,"random_low_H3K27ac"] <- sum(proximal_random_low[match(grouping[[group]],proximal_random_low[,1]),2])
        mat_barplot_chipseq_proximal[group,"random_low_H3K27ac_chipseq"] <- sum(proximal_random_low_chipseq[match(grouping[[group]],proximal_random_low[,1]),2])
}


bp_coords <- barplot(mat_barplot_proximal, beside=FALSE, cex.main=0.9, col=unlist(cols_list), las=2, main=paste(cell, "proximal"), xaxt="n", ylim=c(0,1), cex.axis=0.6,border=NA,space=c(2,2,2,2))
barplot(as.matrix(mat_barplot_chipseq_proximal[,"dnase_high_H3K27ac_chipseq"]),  space=bp_coords[1]+0.5, beside=FALSE, cex.main=0.9, col=unlist(cols_list), axes=FALSE, add=TRUE, border=NA)
barplot(as.matrix(mat_barplot_chipseq_proximal[,"dnase_low_H3K27ac_chipseq"]),  space=bp_coords[2]+0.5, beside=FALSE, cex.main=0.9, col=unlist(cols_list), axes=FALSE, add=TRUE, border=NA)
barplot(as.matrix(mat_barplot_chipseq_proximal[,"random_high_H3K27ac_chipseq"]),  space=bp_coords[3]+0.5, beside=FALSE, cex.main=0.9, col=unlist(cols_list), axes=FALSE, add=TRUE, border=NA)
barplot(as.matrix(mat_barplot_chipseq_proximal[,"random_low_H3K27ac_chipseq"]),  space=bp_coords[4]+0.5, beside=FALSE, cex.main=0.9, col=unlist(cols_list), axes=FALSE, add=TRUE, border=NA)
text(x=bp_coords+0.5,  par("usr")[3]-0.01, labels=colnames(mat_barplot_proximal), cex=0.6, srt = 30, pos = 2, xpd = TRUE)
#legend(15,0.8, c(names(cols_list),"fraction with ChIP-seq intersection"), col=c(unlist(cols_list),"black"), pch=15, bty="n", cex=0.8, pt.cex=1.2)
legend(15,0.8, names(cols_list), col=unlist(cols_list), pch=15, bty="n", cex=0.8, pt.cex=1.2)
dev.off()
