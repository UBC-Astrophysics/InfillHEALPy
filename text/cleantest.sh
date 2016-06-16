ctioga2 --name cleantest cleantest_sorted.dat \
	--xlabel 'Pearson Correlation Coefficient $r$' \
        --ylabel 'Fraction of Trials with the Measured $r_m < r$' \
	--math /xrange 0.05:0.55 'erf((x-0.26715063626373625)/sqrt(2)/0.099593036538962398)*0.5+0.5' 

