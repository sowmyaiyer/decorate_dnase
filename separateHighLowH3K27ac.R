cell <- commandArgs(TRUE)[1]
random_distal_H3K27ac_table <- as.data.frame(read.table(paste("mean_H3K27ac.",cell,".randombg3.distal.txt",sep="")))
dnase_distal_H3K27ac_table <- as.data.frame(read.table(paste("mean_H3K27ac.",cell,".distal.txt",sep="")))
total_bp_in_distal_dnase_peaks <- sum(as.numeric(dnase_distal_H3K27ac_table[,3])-as.numeric(dnase_distal_H3K27ac_table[,2]))
total_bp_in_distal_random_peaks <- sum(as.numeric(random_distal_H3K27ac_table[,3])-as.numeric(random_distal_H3K27ac_table[,2]))

random_proximal_H3K27ac_table <- as.data.frame(read.table(paste("mean_H3K27ac.",cell,".randombg3.proximal.txt",sep="")))
dnase_proximal_H3K27ac_table <- as.data.frame(read.table(paste("mean_H3K27ac.",cell,".proximal.txt",sep="")))
total_bp_in_proximal_dnase_peaks <- sum(as.numeric(dnase_proximal_H3K27ac_table[,3])-as.numeric(dnase_proximal_H3K27ac_table[,2]))
total_bp_in_proximal_random_peaks <- sum(as.numeric(random_proximal_H3K27ac_table[,3])-as.numeric(random_proximal_H3K27ac_table[,2]))

bg_H3K27_distal_values <- as.numeric(random_distal_H3K27ac_table[,ncol(random_distal_H3K27ac_table)])
iles_5 <- quantile(bg_H3K27_distal_values, prob=seq(0,1,0.05))
distal_threshold <- as.numeric(iles_5[length(iles_5)-1]) # last ile

bg_H3K27_proximal_values <- as.numeric(random_proximal_H3K27ac_table[,ncol(random_proximal_H3K27ac_table)])
iles_5 <- quantile(bg_H3K27_proximal_values, prob=seq(0,1,0.05))
proximal_threshold <- as.numeric(iles_5[length(iles_5)-1]) # last ile

write.table(x=dnase_distal_H3K27ac_table[which(as.numeric(dnase_distal_H3K27ac_table[,5]) >= distal_threshold),],file=paste("dnase_",cell,".combined.merged.distal.high_H3K27ac.narrowPeak",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(x=dnase_distal_H3K27ac_table[which(as.numeric(dnase_distal_H3K27ac_table[,5]) < distal_threshold),],file=paste("dnase_",cell,".combined.merged.distal.low_H3K27ac.narrowPeak",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(x=random_distal_H3K27ac_table[which(as.numeric(random_distal_H3K27ac_table[,5]) >= distal_threshold),],file=paste("random_",cell,".combined.merged.distal.high_H3K27ac.narrowPeak",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(x=random_distal_H3K27ac_table[which(as.numeric(random_distal_H3K27ac_table[,5]) < distal_threshold),],file=paste("random_",cell,".combined.merged.distal.low_H3K27ac.narrowPeak",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

write.table(x=dnase_proximal_H3K27ac_table[which(as.numeric(dnase_proximal_H3K27ac_table[,5]) >= proximal_threshold),],file=paste("dnase_",cell,".combined.merged.proximal.high_H3K27ac.narrowPeak",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(x=dnase_proximal_H3K27ac_table[which(as.numeric(dnase_proximal_H3K27ac_table[,5]) < proximal_threshold),],file=paste("dnase_",cell,".combined.merged.proximal.low_H3K27ac.narrowPeak",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(x=random_proximal_H3K27ac_table[which(as.numeric(random_proximal_H3K27ac_table[,5]) >= proximal_threshold),],file=paste("random_",cell,".combined.merged.proximal.high_H3K27ac.narrowPeak",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(x=random_proximal_H3K27ac_table[which(as.numeric(random_proximal_H3K27ac_table[,5]) < proximal_threshold),],file=paste("random_",cell,".combined.merged.proximal.low_H3K27ac.narrowPeak",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

