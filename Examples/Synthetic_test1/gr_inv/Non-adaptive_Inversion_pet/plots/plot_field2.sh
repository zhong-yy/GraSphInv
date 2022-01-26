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

gmt begin synthetic_test_fixed_mesh_gr_misfit jpg,eps E300
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 1p
	gmt set FORMAT_GEO_MAP ddd:mm:ssF
        gmt basemap -JL100/30/25/35/4.2c -R$Region1 -BWSen -Bafg #-Bxya5f5g5
        gmt grd2cpt $OBSERVED_FIELD"_grd" -Chaxby -E20   #-Z
        gmt grdimage  $OBSERVED_FIELD"_grd"  -E2000 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"mGal"
        echo "(a) Synthetic g@-r@- anomalies" | gmt text -F+cTL+f9p, -D-0.9c/0.8c -N
        
        gmt basemap -JL100/30/25/35/4.2c -R$Region1 -BWSen -Bafg -X+6.5c #-Bxya5f5g5
        gmt grd2cpt $PREDICTED_FIELD"_grd" -Chaxby -E20 #-Z
        gmt grdimage  $PREDICTED_FIELD"_grd"  -E1000 -Q 
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"mGal"
        echo "(b) Predicted g@-r@- anomalies" | gmt text -F+cTL+f9p, -D-0.9c/0.8c -N
        
        gmt basemap -JL100/30/25/35/4.2c -R$Region1 -BWSen -Bafg -X-6.5c -Y-7.5c #-Bxya5f5g5
        gmt grd2cpt "residual_grd" -Chaxby -E20
        gmt grdimage  "residual_grd"  -E1000 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxa10f5 -By+L"mGal"
        echo "(c) Residuals" | gmt text -F+cTL+f9p, -D-0.9c/0.85c -N        	        	        
        
		gmt basemap  -JX4.3c/4.6c -R-38/38/0/450 -BWSrt -Bxa20f5+l"Residual value (mGal)" -Bya100f100+l"Counts" -X6.5c -Y-0.2c
		gmt histogram "residual_gmt_format"   -W2 -L1p -i2
#		gmt histogram "residual_gmt_format" -JX? -Bxa+l"mGal" -Bya200+l"Counts" -BWSrt -W1+b -L1p -i2
		echo "(d) Histogram of the residuals" | gmt text -F+cTL+f9p, -D-0.9c/0.87c -N
		meanstd=`python get_mean_std.py residual_gmt_format`
		echo "${meanstd}" | gmt text -F+cTL+f9p, -D0.1c/-0.2c -N 	
gmt end show
rm $OBSERVED_FIELD"_gmt_format"
rm $OBSERVED_FIELD"_grd"
rm $PREDICTED_FIELD"_gmt_format"
rm $PREDICTED_FIELD"_grd"
rm residual_grd
rm residual_gmt_format
