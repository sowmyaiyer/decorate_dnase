type=$1
echo $type
if [[ -f $HOME/nearline/decorate_dnase/alldnase/all_dnase_peaks.bed ]]; then
        mv $HOME/nearline/decorate_dnase/alldnase/all_dnase_peaks.bed ~/farline/RECYCLE_BIN/
fi
if [[ -f $HOME/nearline/decorate_dnase/${type}_excludable.bed ]]; then
	mv $HOME/nearline/decorate_dnase/${type}_excludable.bed  ~/farline/RECYCLE_BIN/
fi

for peakGz in $HOME/nearline/decorate_dnase/alldnase/ENCSR*DNase*narrowPeak.gz
do
        echo $peakGz
	nm=`basename $peakGz | sed 's/\.gz//g'`
        zcat $peakGz > $HOME/nearline/decorate_dnase/alldnase/$nm
	#checking for malformed peak files. may not be needed
#        awk -F"\t" '{ if (NF != 10) print ARGV[1] }' $nm
        awk -F"\t" '{ print $1"\t"$2"\t"$3 }' $HOME/nearline/decorate_dnase/alldnase/$nm >> $HOME/nearline/decorate_dnase/alldnase/all_dnase_peaks.bed
done

#combined file would be very large. Split into chroms, sort each chrom then concat
awk '{ print $0 >> "/home/si14w/nearline/decorate_dnase/alldnase/alldnase."$1".bed" }' $HOME/nearline/decorate_dnase/alldnase/all_dnase_peaks.bed
awk '{ print $1"\t"$2"\t"$3 >> "/home/si14w/nearline/decorate_dnase/alldnase/alldnase."$1".bed"}' $HOME/nearline/decorate_dnase/wgEncodeDacMapabilityConsensusExcludable.bed
## Include the following line if creating exclusion file for distant random bg regions
if [[ $type == "distal" ]]; then
	echo distal
	awk '{ print $1"\t"$2"\t"$3 >> "/home/si14w/nearline/decorate_dnase/alldnase/alldnase."$1".bed"}' $HOME/nearline/enhancer_predictions/gencodeV19_2K_promoter.bed
fi
for f in $HOME/nearline/decorate_dnase/alldnase/alldnase.chr*.bed
do
#       echo """
        sort -k1,1V -k2,2n $f > $f.sorted.${type}_excludable
#       """ > sort.$f.bsub
done

for chrFile in `ls $HOME/nearline/decorate_dnase/alldnase/*sorted.${type}_excludable | awk '{ print $NF}' | sort -k1,1V`
do
        echo $chrFile
        cat $chrFile >> $HOME/nearline/decorate_dnase/${type}_excludable.bed
done

#rm $HOME/nearline/decorate_dnase/alldnase/ENCSR*DNase*narrowPeak
#rm $HOME/nearline/decorate_dnase/alldnase/alldnase.chr*.bed
#rm $HOME/nearline/decorate_dnase/alldnase/all_dnase_peaks.bed
