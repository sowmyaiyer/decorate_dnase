for bigWigFile in ./ENC*bigWig
do
	cellLine=`basename $bigWigFile | awk -F"_" '{ print $5 }'`
#	for dnasepeak in `ls alldnase/*$cellLine*narrowPeak.gz`
#	do
#		zcat $dnasepeak >> dnase_$cellLine.combined.narrowPeak
#	done
#	sort -k1,1 -k2,2n dnase_$cellLine.combined.narrowPeak | bedtools merge -i stdin > dnase_$cellLine.combined.merged.narrowPeak
echo """	
	bedtools intersect -a dnase_$cellLine.combined.merged.narrowPeak -b /home/si14w/nearline/enhancer_predictions/gencodeV19_2K_promoter.bed  -u > dnase_$cellLine.combined.merged.proximal.narrowPeak
	bedtools intersect -a dnase_$cellLine.combined.merged.narrowPeak -b /home/si14w/nearline/enhancer_predictions/gencodeV19_2K_promoter.bed  -v > dnase_$cellLine.combined.merged.distal.narrowPeak

	sort -k1,1 -k2,2n dnase_$cellLine.combined.merged.distal.narrowPeak | awk '{ print \$0\"\\t\"NR}' | bigWigAverageOverBed $bigWigFile stdin out.$cellLine.tab -bedOut=mean_H3K27ac.$cellLine.distal.txt
	sort -k1,1 -k2,2n dnase_$cellLine.combined.merged.proximal.narrowPeak | awk '{ print \$0\"\\t\"NR}'  | bigWigAverageOverBed $bigWigFile stdin out.$cellLine.tab -bedOut=mean_H3K27ac.$cellLine.proximal.txt

	bedtools  shuffle -i dnase_$cellLine.combined.merged.distal.narrowPeak -g $HOME/TR/hg19.genome -excl distal_excludable.bed  > dnase_$cellLine.combined.merged.distal.random3.narrowPeak

        bedtools  shuffle -i dnase_$cellLine.combined.merged.proximal.narrowPeak -g $HOME/TR/hg19.genome -incl /home/si14w/nearline/enhancer_predictions/gencodeV19_2K_promoter.bed -excl proximal_excludable.bed > dnase_$cellLine.combined.merged.proximal.random3.narrowPeak

        sort -k1,1 -k2,2n dnase_$cellLine.combined.merged.distal.random3.narrowPeak | awk '{ print \$0\"\\t\"NR}' | bigWigAverageOverBed $bigWigFile stdin out2.$cellLine.tab -bedOut=mean_H3K27ac.$cellLine.randombg3.distal.txt
        sort -k1,1 -k2,2n dnase_$cellLine.combined.merged.proximal.random3.narrowPeak | awk '{ print \$0\"\\t\"NR}' | bigWigAverageOverBed $bigWigFile stdin out2.$cellLine.tab -bedOut=mean_H3K27ac.$cellLine.randombg3.proximal.txt

	Rscript drawH3K27Ac.R ${cellLine}
""" > H3K27Ac.${cellLine}.bsub
done
