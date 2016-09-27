#       --draw-line 0.0205115,0 0.0205115,1 \
#	--math /xrange 0.005:0.035 \
#	'(erf( (x-0.0140924)/0.00270706/sqrt(2))+1)/2'  \
ctioga2 --name dhzoa \
	--setup-grid 1x2 --inset grid:0,0 \
	dist_hzoa.txt  \
	--math /xrange 0.018:0.032 \
	--draw-line 0.0284629,0 0.0284629,1 \
	'(erf( (x-0.0242999)/0.00215178/sqrt(2))+1)/2'  \
        --xlabel 'Mean of Map on Random Points' \
        --ylabel '$C( < m)$' 
