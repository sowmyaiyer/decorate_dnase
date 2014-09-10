for cell in {"GM12878","HepG2","K562","endothelial.cell.of.umbilical.vein","H1-hESC","HeLa-S3"}
do
	echo $cell
	bedtools intersect -a dnase_${cell}.combined.merged.distal.high_H3K27ac.narrowPeak -b allchipseq/all_chipseq_peaks.${cell}.bed.sorted.merged -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9 -o min | awk '{ print $NF}' > ${cell}_dnase_high_H3K27ac_intersect_mergedChipSeq.distal.txt
	bedtools intersect -a dnase_${cell}.combined.merged.distal.low_H3K27ac.narrowPeak -b allchipseq/all_chipseq_peaks.${cell}.bed.sorted.merged -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9 -o min  | awk '{ print $NF}' > ${cell}_dnase_low_H3K27ac_intersect_mergedChipSeq.distal.txt	
	bedtools intersect -a random_${cell}.combined.merged.distal.high_H3K27ac.narrowPeak -b allchipseq/all_chipseq_peaks.${cell}.bed.sorted.merged -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9 -o min  | awk '{ print $NF}'> ${cell}_random_high_H3K27ac_intersect_mergedChipSeq.distal.txt
	bedtools intersect -a random_${cell}.combined.merged.distal.low_H3K27ac.narrowPeak -b allchipseq/all_chipseq_peaks.${cell}.bed.sorted.merged -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9 -o min  | awk '{ print $NF}'> ${cell}_random_low_H3K27ac_intersect_mergedChipSeq.distal.txt	
	
	bedtools intersect -a dnase_${cell}.combined.merged.proximal.high_H3K27ac.narrowPeak -b allchipseq/all_chipseq_peaks.${cell}.bed.sorted.merged -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9 -o min  | awk '{ print $NF}' > ${cell}_dnase_high_H3K27ac_intersect_mergedChipSeq.proximal.txt
 	
	bedtools intersect -a dnase_${cell}.combined.merged.proximal.low_H3K27ac.narrowPeak -b allchipseq/all_chipseq_peaks.${cell}.bed.sorted.merged -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9 -o min  | awk '{ print $NF}' > ${cell}_dnase_low_H3K27ac_intersect_mergedChipSeq.proximal.txt
	bedtools intersect -a random_${cell}.combined.merged.proximal.high_H3K27ac.narrowPeak -b allchipseq/all_chipseq_peaks.${cell}.bed.sorted.merged -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9 -o min  | awk '{ print $NF}' > ${cell}_random_high_H3K27ac_intersect_mergedChipSeq.proximal.txt

	bedtools intersect -a random_${cell}.combined.merged.proximal.low_H3K27ac.narrowPeak -b allchipseq/all_chipseq_peaks.${cell}.bed.sorted.merged -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9 -o min  | awk '{ print $NF}' > ${cell}_random_low_H3K27ac_intersect_mergedChipSeq.proximal.txt 	
	
	Rscript chromHmmIntersection_high_low_K27ac_version2.R ${cell}
done
