#	hzoa_cum.dat --yscale 0.00113481616 @1:3 --draw-line 0.0,0.0 1,1 \
ctioga2 --name hzoa_cum \
	test_dist.dat \
	--xscale 0.00004562043796 --yscale 0.004310344828 \
	hzoa_cum.dat \
	--draw-line 0,0 1,1 \
	--draw-line 0,0.5 1,0.5 \
	--draw-line 0,0.8 1,0.8 \
	--draw-line 0.47,0 0.47,1 \
	--draw-line 0,0.8 1,0.8 \
	--draw-line 0.725,0 0.725,1 \
	-x 'Fraction of Solid Angle Surveyed' -y 'Fraction of Galaxies Found'
