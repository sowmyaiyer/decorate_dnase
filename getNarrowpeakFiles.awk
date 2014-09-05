BEGIN {FS="\t"}
{
	if ($7==cellLine)
        {
		wget=sprintf("wget --background %s -O $HOME/nearline/decorate_dnase/%s.narrowPeak.gz",$14,$24)
		print(wget)
	}
}
