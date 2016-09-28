# @'$1:($3*0.25+$1*0.75)' \
ctioga2 --name test_dist test_dist.dat @1:3 \
        --draw-line 0,0 1,1 \
	-x 'Fraction of Solid Angle Surveyed' \
	-y 'Fraction of Galaxies Found'

