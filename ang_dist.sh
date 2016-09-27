awk '($1>10 || $1<-10)' ang_res.dat | sort -g -k 2 | awk '{print $2, NR/3399.,NR,$0}' > ang_data_sorted.dat
ctioga2 --name ang_dist ang_data_sorted.dat \
        --xlabel 'Pearson Correlation Coefficient $r$' \
        --ylabel 'Fraction of Trials with the Measured $r_m < r$' \
	--math /xrange -0.2:0.2 '(erf((x+0.010293668411297328)/0.043220082133824014/sqrt(2))+1)/2.0' \
