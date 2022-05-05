Trr=../dobs_T_rr
Trtheta=../dobs_T_rtheta
Trphi=../dobs_T_rphi
Tthetatheta=../dobs_T_thetatheta
Tthetaphi=../dobs_T_thetaphi

Trr_predicted=../dpredicted_T_rr
Trtheta_predicted=../dpredicted_T_rtheta
Trphi_predicted=../dpredicted_T_rphi
Tthetatheta_predicted=../dpredicted_T_thetatheta
Tthetaphi_predicted=../dpredicted_T_thetaphi

python swap_2cols_of_Trr_data.py $Trr $Trr_predicted "_gmt_format"
python swap_2cols_of_Trr_data.py $Trtheta $Trtheta_predicted "_gmt_format"
python swap_2cols_of_Trr_data.py $Trphi $Trphi_predicted "_gmt_format"
python swap_2cols_of_Trr_data.py $Tthetatheta $Tthetatheta_predicted "_gmt_format"
python swap_2cols_of_Trr_data.py $Tthetaphi $Tthetaphi_predicted "_gmt_format"

Region1=90/110/20/40
xinc_yinc=0.5/0.5


gmt xyz2grd $Trr"_gmt_format" -I$xinc_yinc -R$Region1 -G$Trr"_grd"
gmt xyz2grd $Trtheta"_gmt_format" -I$xinc_yinc -R$Region1 -G$Trtheta"_grd"
gmt xyz2grd $Trphi"_gmt_format" -I$xinc_yinc -R$Region1 -G$Trphi"_grd"
gmt xyz2grd $Tthetatheta"_gmt_format" -I$xinc_yinc -R$Region1 -G$Tthetatheta"_grd"
gmt xyz2grd $Tthetaphi"_gmt_format" -I$xinc_yinc -R$Region1 -G$Tthetaphi"_grd"


gmt xyz2grd $Trr_predicted"_gmt_format" -I$xinc_yinc -R$Region1 -G$Trr_predicted"_grd"
gmt xyz2grd $Trtheta_predicted"_gmt_format" -I$xinc_yinc -R$Region1 -G$Trtheta_predicted"_grd"
gmt xyz2grd $Trphi_predicted"_gmt_format" -I$xinc_yinc -R$Region1 -G$Trphi_predicted"_grd"
gmt xyz2grd $Tthetatheta_predicted"_gmt_format" -I$xinc_yinc -R$Region1 -G$Tthetatheta_predicted"_grd"
gmt xyz2grd $Tthetaphi_predicted"_gmt_format" -I$xinc_yinc -R$Region1 -G$Tthetaphi_predicted"_grd"


gmt xyz2grd $Trr"_residual_gmt_format" -I$xinc_yinc -R$Region1 -G$Trr"_residual_grd"
gmt xyz2grd $Trtheta"_residual_gmt_format" -I$xinc_yinc -R$Region1 -G$Trtheta"_residual_grd"
gmt xyz2grd $Trphi"_residual_gmt_format" -I$xinc_yinc -R$Region1 -G$Trphi"_residual_grd"
gmt xyz2grd $Tthetatheta"_residual_gmt_format" -I$xinc_yinc -R$Region1 -G$Tthetatheta"_residual_grd"
gmt xyz2grd $Tthetaphi"_residual_gmt_format" -I$xinc_yinc -R$Region1 -G$Tthetaphi"_residual_grd"


#lon_inv=(67.5 109.5)
#lat_inv=(21.5 47.5)

#inversion_region='inversion_region'

#python make_frame.py ${lon_inv[0]} ${lon_inv[1]} ${lat_inv[0]} ${lat_inv[1]} $inversion_region
gmt begin ggt_observed pdf E300
	gmt set FONT_TAG 9p
	gmt set FONT_TITLE 9p
	gmt set MAP_TITLE_OFFSET 9p 
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 1p,black
	gmt set FORMAT_GEO_MAP ddd:mm:ssF
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg #-Bxya5f5g5
        gmt grd2cpt $Trr"_grd" -Chaxby -E20
        gmt grdimage  $Trr"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(a) Synthetic "T@-rr@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
        
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c #-Bxya5f5g5
        gmt grd2cpt $Trtheta"_grd" -Chaxby -E20
        gmt grdimage  $Trtheta"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(b) Synthetic "T@-r@~'\161'@~@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
        
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c #-Bxya5f5g5
        gmt grd2cpt $Trphi"_grd" -Chaxby -E20
        gmt grdimage  $Trphi"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(c) Synthetic "T@-r@~'\152'@~@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
            
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X-5.6c -Y-7.3c #-Bxya5f5g5
        gmt grd2cpt $Tthetatheta"_grd" -Chaxby -E20
        gmt grdimage  $Tthetatheta"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(d) Synthetic "T@-@~'\161''\161'@~@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N        	        	        
        
		gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c #-Bxya5f5g5
		gmt grd2cpt $Tthetaphi"_grd" -Chaxby -E20
        gmt grdimage  $Tthetaphi"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(d) Synthetic "T@-@~'\161''\152'@~@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N    
gmt end show

gmt begin ggt_predicted pdf E300
	gmt set FONT_TAG 9p
	gmt set FONT_TITLE 9p
	gmt set MAP_TITLE_OFFSET 9p 
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 1p,black
	gmt set FORMAT_GEO_MAP ddd:mm:ssF
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg #-Bxya5f5g5
        gmt grd2cpt $Trr_predicted"_grd" -Chaxby -E20
        gmt grdimage  $Trr_predicted"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(a) Predicted "T@-rr@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
        
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c #-Bxya5f5g5
        gmt grd2cpt $Trtheta_predicted"_grd" -Chaxby -E20
        gmt grdimage  $Trtheta_predicted"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(b) Predicted "T@-r@~'\161'@~@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
        
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c #-Bxya5f5g5
        gmt grd2cpt $Trphi_predicted"_grd" -Chaxby -E20
        gmt grdimage  $Trphi_predicted"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(c) Predicted "T@-r@~'\152'@~@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
            
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X-5.6c -Y-7.3c #-Bxya5f5g5
        gmt grd2cpt $Tthetatheta_predicted"_grd" -Chaxby -E20
        gmt grdimage  $Tthetatheta_predicted"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(d) Predicted "T@-@~'\161''\161'@~@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N        	        	        
        
		gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c #-Bxya5f5g5
		gmt grd2cpt $Tthetaphi_predicted"_grd" -Chaxby -E20
        gmt grdimage  $Tthetaphi_predicted"_grd"  -E300 -Q
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(e) Predicted "T@-@~'\161''\152'@~@-" anomalies" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N    
gmt end show

gmt begin ggt_residuals pdf E300
	gmt set FONT_TAG 9p
	gmt set FONT_TITLE 9p
	gmt set MAP_TITLE_OFFSET 9p 
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 1p,black
	gmt set FORMAT_GEO_MAP ddd:mm:ssF
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg
        gmt grd2cpt  $Trr"_residual_grd" -Chaxby -E20
        gmt grdimage  $Trr"_residual_grd"  -E300 -Q -Chaxby
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(a) Residuals of  "T@-rr@-" data" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
        
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c
        gmt grd2cpt  $Trtheta"_residual_grd" -Chaxby -E20
        gmt grdimage  $Trtheta"_residual_grd"  -E300 -Q -Chaxby
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(b) Residuals of "T@-r@~'\161'@~@-" data" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
        
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c
        gmt grd2cpt  $Trphi"_residual_grd" -Chaxby -E20
        gmt grdimage  $Trphi"_residual_grd"  -E300 -Q -Chaxby
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(c) Residuals of "T@-r@~'\152'@~@-" data" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N
            
        gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X-5.6c -Y-7.3c
        gmt grd2cpt  $Tthetatheta"_residual_grd" -Chaxby -E20
        gmt grdimage  $Tthetatheta"_residual_grd"  -E300 -Q -Chaxby
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(d) Residuals of "T@-@~'\161''\161'@~@-" data" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N        	        	        
        
		gmt basemap -JL100/30/25/35/3.75c -R$Region1 -BWSen -Bafg -X+5.6c
		gmt grd2cpt  $Tthetaphi"_residual_grd" -Chaxby -E20
        gmt grdimage  $Tthetaphi"_residual_grd"  -E300 -Q -Chaxby
        gmt colorbar -DJBC+h+w3.8c/0.2c+o0.c/0.9c -Bxaf -By+L"E"
        echo "(e) Residuals of "T@-@~'\161''\152'@~@-" data" | gmt text -F+cTL+f9p, -D-0.8c/0.85c -N    
gmt end show

gmt begin ggt_residuals_histogram pdf E300
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 1p,black
	gmt set FORMAT_GEO_MAP ddd:mm:ssF
        #gmt basemap -JL100/30/25/35/4.2c -R$Region1 -BWSen -Bafg #-Bxya5f5g5
        gmt basemap  -JX4.25c/3.8c -R-1.8/1.8/0/500 -BWSrt -Bxa1f0.2+l"Residual value (E)" -Bya200f100+l"Counts"
        gmt histogram $Trr"_residual_gmt_format"   -W0.1 -L1p -i2
		meanstd=`python get_mean_std.py "$Trr"_residual_gmt_format`
        echo "(a) Residuals of  "T@-rr@-" data" | gmt text -F+cTL+f9p, -D-0.6c/0.6c -N
        echo "${meanstd}" | gmt text -F+cTL+f9p, -D0.05c/-0.2c -N
        
        gmt basemap  -JX4.25c/3.8c -R-1.8/1.8/0/500 -BWSrt -Bxa1f0.2+l"Residual value (E)" -Bya200f100+l"Counts" -X+6.5c
        gmt histogram $Trtheta"_residual_gmt_format"   -W0.1 -L1p -i2
		meanstd=`python get_mean_std.py "$Trtheta"_residual_gmt_format`
        echo "(b) Residuals of "T@-r@~'\161'@~@-" data" | gmt text -F+cTL+f9p, -D-0.6c/0.6c -N
        echo "${meanstd}" | gmt text -F+cTL+f9p, -D0.05c/-0.2c -N
        
        gmt basemap  -JX4.25c/3.8c -R-1.8/1.8/0/500 -BWSrt -Bxa1f0.2+l"Residual value (E)" -Bya200f100+l"Counts" -X+6.5c
        gmt histogram $Trphi"_residual_gmt_format"   -W0.1 -L1p -i2
		meanstd=`python get_mean_std.py "$Trphi"_residual_gmt_format`        
        echo "(c) Residuals of "T@-r@~'\152'@~@-" data" | gmt text -F+cTL+f9p, -D-0.6c/0.6c -N    
        echo "${meanstd}" | gmt text -F+cTL+f9p, -D0.05c/-0.2c -N  
        
        gmt basemap  -JX4.25c/3.8c -R-1.8/1.8/0/500 -BWSrt -Bxa1f0.2+l"Residual value (E)" -Bya200f100+l"Counts" -X-6.5c -Y-6.5c
        gmt histogram $Tthetatheta"_residual_gmt_format"   -W0.1 -L1p -i2
		meanstd=`python get_mean_std.py "$Tthetatheta"_residual_gmt_format`
        echo "(d) Residuals of "T@-@~'\161''\161'@~@-" data" | gmt text -F+cTL+f9p, -D-0.6c/0.6c -N
        echo "${meanstd}" | gmt text -F+cTL+f9p, -D0.05c/-0.2c -N        
  	        	        
        
        gmt basemap  -JX4.25c/3.8c -R-1.8/1.8/0/500 -BWSrt -Bxa1f0.2+l"Residual value (E)" -Bya200f100+l"Counts" -X+6.5c
        gmt histogram $Tthetaphi"_residual_gmt_format"   -W0.1 -L1p -i2
		meanstd=`python get_mean_std.py "$Tthetaphi"_residual_gmt_format`
        echo "(e) Residuals of "T@-@~'\161''\152'@~@-" data" | gmt text -F+cTL+f9p, -D-0.6c/0.6c -N
        echo "${meanstd}" | gmt text -F+cTL+f9p, -D0.05c/-0.2c -N
gmt end show

rm $Trr"_gmt_format"
rm $Trr"_grd"
rm $Trr_predicted"_gmt_format"
rm $Trr_predicted"_grd"
rm $Trr"_residual_gmt_format"
rm $Trr"_residual_grd"

rm $Trtheta"_gmt_format"
rm $Trtheta"_grd"
rm $Trtheta_predicted"_gmt_format"
rm $Trtheta_predicted"_grd"
rm $Trtheta"_residual_gmt_format"
rm $Trtheta"_residual_grd"

rm $Tthetatheta"_gmt_format"
rm $Tthetatheta"_grd"
rm $Tthetatheta_predicted"_gmt_format"
rm $Tthetatheta_predicted"_grd"
rm $Tthetatheta"_residual_gmt_format"
rm $Tthetatheta"_residual_grd"

rm $Trphi"_gmt_format"
rm $Trphi"_grd"
rm $Trphi_predicted"_gmt_format"
rm $Trphi_predicted"_grd"
rm $Trphi"_residual_gmt_format"
rm $Trphi"_residual_grd"

rm $Tthetaphi"_gmt_format"
rm $Tthetaphi"_grd"
rm $Tthetaphi_predicted"_gmt_format"
rm $Tthetaphi_predicted"_grd"
rm $Tthetaphi"_residual_gmt_format"
rm $Tthetaphi"_residual_grd"
