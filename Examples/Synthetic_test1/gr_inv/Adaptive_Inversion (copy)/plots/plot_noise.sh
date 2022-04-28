OBSERVED_FIELD=../dobs_g_r
EXACT_FIELD=./dobs_exact_g_r

echo $OBSERVED_FIELD
echo $EXACT_FIELD

python swap_2cols_of_g_data.py $OBSERVED_FIELD $EXACT_FIELD "_gmt_format"

Region1=90/110/20/40
xinc_yinc=0.5/0.5
gmt xyz2grd $OBSERVED_FIELD"_gmt_format" -I$xinc_yinc -R$Region1 -G$OBSERVED_FIELD"_grd"
gmt xyz2grd $EXACT_FIELD"_gmt_format" -I$xinc_yinc -R$Region1 -G$EXACT_FIELD"_grd"
gmt xyz2grd "residual_gmt_format" -I$xinc_yinc -R$Region1 -G"residual_grd"


#lon_inv=(67.5 109.5)
#lat_inv=(21.5 47.5)

#inversion_region='inversion_region'

#python make_frame.py ${lon_inv[0]} ${lon_inv[1]} ${lat_inv[0]} ${lat_inv[1]} $inversion_region

gmt begin noise jpg,eps E300
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set FORMAT_GEO_MAP ddd:mm:ssF        
        gmt basemap -JL100/30/25/35/4.2c -R$Region1 -BWSen -Bafg
        #gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5	        
        #gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
        gmt grd2cpt "residual_grd" -Chaxby -E20
        gmt grdimage  "residual_grd"  -E1000 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"mGal"
        echo "(a) Residuals" | gmt text -F+cTL+f9p, -D-0.9c/0.85c -N        	        	        
        
		gmt basemap  -JX4.3c/4.8c -R-35/35/0/500 -BWSrt -Bxa20f5+l"Residual value (mGal)" -Bya100f100+l"Counts" -X6.3c -Y-0.2c
		gmt histogram "residual_gmt_format"   -W2 -L1p -i2
#		gmt histogram "residual_gmt_format" -JX? -Bxa+l"mGal" -Bya200+l"Counts" -BWSrt -W1+b -L1p -i2
		echo "(b) Histogram of the residuals" | gmt text -F+cTL+f9p, -D-0.73c/0.72c -N
		meanstd=`python get_mean_std.py residual_gmt_format`
		echo "${meanstd}" | gmt text -F+cTL+f9p, -D0.1c/-0.2c -N 	
	        #gmt basemap -JL88/34/30/43/? -R69/107/23/45 -BWSen -Bxya5f5g5
	        #gmt basemap -JM88/34/? -R69/107/23/45 -BwSEn -Baf5g5	        
	        #gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
	        #gmt grdimage  $ISOSTASY_EFFECT"_grd"  -E200 -t30 -Q -Cpanoply
	        #gmt colorbar -DJBC+h+w7c/0.2c+o0.c/1c -Bxa100f -By+L"mGal"
gmt end show
rm $OBSERVED_FIELD"_gmt_format"
rm $OBSERVED_FIELD"_grd"
rm $EXACT_FIELD"_gmt_format"
rm $EXACT_FIELD"_grd"
