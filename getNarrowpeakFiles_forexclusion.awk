BEGIN {FS="\t"}
{
	if ($11 == "bed_narrowPeak")
	{
		wget=sprintf("wget --background %s -O $HOME/nearline/decorate_dnase/alldnase/%s.narrowPeak.gz",$14,$24)
		print(wget)
	}
}
