for cell in {"GM12878","HepG2","K562","endothelial.cell.of.umbilical.vein","H1-hESC","HeLa-S3"}
do
	echo $cell
	Rscript getBackgroundStats.R ${cell}
	
	bedtools intersect -a dnase_${cell}.combined.merged.distal.high_H3K27ac.narrowPeak -b wgEncodeAwgSegmentationChromhmm.${cell}.bed -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9,15 -o collapse,collapse > ${cell}_dnase_high_H3K27ac_intersect_chromHmm.25.states.distal_v2.txt
	awk -f pickState.awk ${cell}_dnase_high_H3K27ac_intersect_chromHmm.25.states.distal_v2.txt | sort |  uniq -c | awk '{ print $2"\t"$1}' > ${cell}_dnase_high_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt

	bedtools intersect -a dnase_${cell}.combined.merged.distal.low_H3K27ac.narrowPeak -b wgEncodeAwgSegmentationChromhmm.${cell}.bed -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9,15 -o collapse,collapse > ${cell}_dnase_low_H3K27ac_intersect_chromHmm.25.states.distal_v2.txt
	awk -f pickState.awk ${cell}_dnase_low_H3K27ac_intersect_chromHmm.25.states.distal_v2.txt | sort |  uniq -c | awk '{ print $2"\t"$1}' > ${cell}_dnase_low_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt

	bedtools intersect -a random_${cell}.combined.merged.distal.high_H3K27ac.narrowPeak -b wgEncodeAwgSegmentationChromhmm.${cell}.bed -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9,15 -o collapse,collapse > ${cell}_random_high_H3K27ac_intersect_chromHmm.25.states.distal_v2.txt
	awk -f pickState.awk ${cell}_random_high_H3K27ac_intersect_chromHmm.25.states.distal_v2.txt | sort |  uniq -c | awk '{ print $2"\t"$1}' > ${cell}_random_high_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt
 	
	bedtools intersect -a random_${cell}.combined.merged.distal.low_H3K27ac.narrowPeak -b wgEncodeAwgSegmentationChromhmm.${cell}.bed -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9,15 -o collapse,collapse > ${cell}_random_low_H3K27ac_intersect_chromHmm.25.states.distal_v2.txt
	awk -f pickState.awk ${cell}_random_low_H3K27ac_intersect_chromHmm.25.states.distal_v2.txt | sort |  uniq -c | awk '{ print $2"\t"$1}' > ${cell}_random_low_H3K27ac_intersect_chromHmm.25.states.distal.counts.txt



	bedtools intersect -a dnase_${cell}.combined.merged.proximal.high_H3K27ac.narrowPeak -b wgEncodeAwgSegmentationChromhmm.${cell}.bed -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9,15 -o collapse,collapse > ${cell}_dnase_high_H3K27ac_intersect_chromHmm.25.states.proximal_v2.txt
	awk -f pickState.awk ${cell}_dnase_high_H3K27ac_intersect_chromHmm.25.states.proximal_v2.txt | sort |  uniq -c | awk '{ print $2"\t"$1}' > ${cell}_dnase_high_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt

	bedtools intersect -a dnase_${cell}.combined.merged.proximal.low_H3K27ac.narrowPeak -b wgEncodeAwgSegmentationChromhmm.${cell}.bed -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9,15 -o collapse,collapse > ${cell}_dnase_low_H3K27ac_intersect_chromHmm.25.states.proximal_v2.txt
	awk -f pickState.awk ${cell}_dnase_low_H3K27ac_intersect_chromHmm.25.states.proximal_v2.txt | sort |  uniq -c | awk '{ print $2"\t"$1}' > ${cell}_dnase_low_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt

	bedtools intersect -a random_${cell}.combined.merged.proximal.high_H3K27ac.narrowPeak -b wgEncodeAwgSegmentationChromhmm.${cell}.bed -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9,15 -o collapse,collapse > ${cell}_random_high_H3K27ac_intersect_chromHmm.25.states.proximal_v2.txt
	awk -f pickState.awk ${cell}_random_high_H3K27ac_intersect_chromHmm.25.states.proximal_v2.txt | sort |  uniq -c | awk '{ print $2"\t"$1}' > ${cell}_random_high_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt
 	
	bedtools intersect -a random_${cell}.combined.merged.proximal.low_H3K27ac.narrowPeak -b wgEncodeAwgSegmentationChromhmm.${cell}.bed -wao | sort -k1,1V -k2,2n | bedtools groupby -c 9,15 -o collapse,collapse > ${cell}_random_low_H3K27ac_intersect_chromHmm.25.states.proximal_v2.txt
	awk -f pickState.awk ${cell}_random_low_H3K27ac_intersect_chromHmm.25.states.proximal_v2.txt | sort |  uniq -c | awk '{ print $2"\t"$1}' > ${cell}_random_low_H3K27ac_intersect_chromHmm.25.states.proximal.counts.txt

	
	for f in ./${cell}_*H3K27ac_intersect_chromHmm.25.states*counts.txt
	do
		join states.26.txt $f -a 1 | awk -F" " '{ if (NF == 1) {print $1"\t"0} else {print $1"\t"$2}}' > $f.tmp
		mv $f.tmp $f
		wc -l $f
	done
	Rscript chromHmmIntersection_high_low_K27ac_version2.R ${cell}
done
