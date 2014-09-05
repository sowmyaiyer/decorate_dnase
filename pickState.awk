{
	if ($4 == ".") 
	{
		print($4)
	} else {
		split($4, states, ",")
		split($5, counts, ",")
		argmax=-1
		max=0
		for (i = 0 ; i< length(counts); i ++)
		{
			if (counts[i] > max)
			{
				max=counts[i]
				argmax=i
			}
		}
		print(states[argmax])
	}
}
