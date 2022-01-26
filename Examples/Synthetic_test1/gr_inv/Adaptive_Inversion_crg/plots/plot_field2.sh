OBSERVED_FIELD=../dobs_g_r
PREDICTED_FIELD=../dpredicted_g_r

echo $OBSERVED_FIELD
echo $PREDICTED_FIELD

python swap_2cols_of_g_data.py $OBSERVED_FIELD $PREDICTED_FIELD "_gmt_format"

Region1=90/110/20/40
xinc_yinc=0.5/0.5
gmt xyz2grd $OBSERVED_FIELD"_gmt_format" -I$xinc_yinc -R$Region1 -G$OBSERVED_FIELD"_grd"
gmt xyz2grd $PREDICTED_FIELD"_gmt_format" -I$xinc_yinc -R$Region1 -G$PREDICTED_FIELD"_grd"
gmt xyz2grd "residual_gmt_format" -I$xinc_yinc -R$Region1 -G"residual_grd"


#lon_inv=(67.5 109.5)
#lat_inv=(21.5 47.5)

#inversion_region='inversion_region'

#python make_frame.py ${lon_inv[0]} ${lon_inv[1]} ${lat_inv[0]} ${lat_inv[1]} $inversion_region

gmt begin synthetic_test_gr_crg_misfit eps E300
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set FORMAT_GEO_MAP ddd:mm:ssF
        gmt basemap -JL100/30/25/35/4.2c -R$Region1 -BWSen -Bafg #-Bxya5f5g5
        gmt grd2cpt $OBSERVED_FIELD"_grd" -Chaxby -E20   #-Z
        gmt grdimage  $OBSERVED_FIELD"_grd"  -E2000 -Q
        gmt colorbar -DJBC+h+w3.9c/0.2c+o0.c/0.9c -Bxaf -By+L"mGal"
        echo "(a) Synthetic g@-r@- anomalies" | gmt text -F+cTL+f10p, -D-0.9c/0.85c -N
        
        gmt basemap -JL100/30/25/35/4.2c -R$Region1 -BWSen -Bafg -X+6.8c #-Bxya5f5g5
        gmt grd2cpt $PREDICTED_FIELD"_grd" -Chaxby -E20 #-Z
        gmt grdimage  $PREDICTED_FIELD"_grd"  -E1000 -Q 
        gmt colorbar -DJBC+h+w3.9c/0.2c+o0.c/0.9c -Bxaf -By+L"mGal"
        echo "(b) Predicted g@-r@- anomalies" | gmt text -F+cTL+f10p, -D-0.9c/0.85c -N
        
        gmt basemap -JL100/30/25/35/4.2c -R$Region1 -BWSen -Bafg -X-6.8c -Y-7.5c #-Bxya5f5g5
        #gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5	        
        #gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
        gmt grd2cpt "residual_grd" -Chaxby -E20
        gmt grdimage  "residual_grd"  -E1000 -Q
        gmt colorbar -DJBC+h+w3.9c/0.2c+o0.c/0.9c -Bxa10f5 -By+L"mGal"
        echo "(c) Residuals" | gmt text -F+cTL+f10p, -D-0.9c/0.85c -N        	        	        
        
		gmt basemap  -JX4.8c/4.6c -R-35/35/0/500 -BWSrt -Bxa20f5+l"Residual value (mGal)" -Bya200f100+l"Counts" -X7.c -Y-0.2c
		gmt histogram "residual_gmt_format"   -W2 -L1p -i2
#		gmt histogram "residual_gmt_format" -JX? -Bxa+l"mGal" -Bya200+l"Counts" -BWSrt -W1+b -L1p -i2
		echo "(d) Histogram of the residuals" | gmt text -F+cTL+f10p, -D-1.08c/0.92c -N
		meanstd=`python get_mean_std.py residual_gmt_format`
		echo "${meanstd}" | gmt text -F+cTL+f10p, -D0.1c/-0.2c -N 	
	        #gmt basemap -JL88/34/30/43/? -R69/107/23/45 -BWSen -Bxya5f5g5
	        #gmt basemap -JM88/34/? -R69/107/23/45 -BwSEn -Baf5g5	        
	        #gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
	        #gmt grdimage  $ISOSTASY_EFFECT"_grd"  -E200 -t30 -Q -Cpanoply
	        #gmt colorbar -DJBC+h+w7c/0.2c+o0.c/1c -Bxa100f -By+L"mGal"
gmt end show
rm $OBSERVED_FIELD"_gmt_format"
rm $OBSERVED_FIELD"_grd"
rm $PREDICTED_FIELD"_gmt_format"
rm $PREDICTED_FIELD"_grd"
rm residual_grd
rm residual_gmt_format
