ctioga2 --yfact 0.001 --name cldist sorted2.lst \
	--math /xrange -0.2:0.2 '(erf(x/0.06624878354105386/sqrt(2))+1)/2.0*1000' \
	--xlabel 'Pearson Correlation Coefficient $r$' \
        --ylabel 'Fraction of Trials with the Measured $r_m < r$' \
