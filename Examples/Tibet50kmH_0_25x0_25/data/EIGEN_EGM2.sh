EIGEN6C4_FIELD=../mantle_g_r
EGM_FIELD=./mantle_g_r_EGM2008

echo $EIGEN6C4_FIELD
echo $EGM_FIELD

python swap_2cols_of_g_data.py $EIGEN6C4_FIELD $EGM_FIELD "_gmt_format"

Region1=65/111/20/48
xinc_yinc=0.5/0.5
gmt xyz2grd $EIGEN6C4_FIELD"_gmt_format" -I$xinc_yinc -R$Region1 -G$EIGEN6C4_FIELD"_grd" -D+xLongitude[degrees_east]+yLatitude[degrees_north]
gmt xyz2grd $EGM_FIELD"_gmt_format" -I$xinc_yinc -R$Region1 -G$EGM_FIELD"_grd" -D+xLongitude[degrees_east]+yLatitude[degrees_north]
gmt xyz2grd "residual_gmt_format" -I$xinc_yinc -R$Region1 -G"residual_grd" -D+xLongitude[degrees_east]+yLatitude[degrees_north]

#lon_inv=(67.5 109.5)
#lat_inv=(21.5 47.5)

#inversion_region='inversion_region'

#python make_frame.py ${lon_inv[0]} ${lon_inv[1]} ${lat_inv[0]} ${lat_inv[1]} $inversion_region

function plot_blocks() {
        # ============绘制板块边界
        #gmt plot CN-plate-neighbor.dat -W1.0p,black -Sf0.5+t+l -G2/138/210
        gmt plot Eurasia_India_PB.dat -W0.9p,black,-
        # ============绘制地块边界
        gmt plot China_tectonic.dat -W0.3p,black,-
        # ============腾冲火山
        #echo 98 25 |gmt plot -St0.25 -Gwhite -W0.3p,black
}

gmt begin comparison_of_EIGEN6C4_EGM2008 eps E300
gmt set FONT 9p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_FRAME_PEN 1.5p,black
gmt set FORMAT_GEO_MAP ddd:mm:ssF
gmt basemap -JQ88/34/5c -R$Region1 -BWSen -Ba10f2g #-Bxya5f5g5
gmt grd2cpt $EIGEN6C4_FIELD"_grd" $EGM_FIELD"_grd" "residual_grd" -Chaxby -E300 -Z
gmt grdimage $EIGEN6C4_FIELD"_grd" -I+d -E300 -Q #-Chaxby
#gmt colorbar -DJBC+h+w3.6c/0.2c+o0.c/0.9c -Bxaf -By+L"mGal"
echo "(a) EIGEN-6C4" | gmt text -F+cTL+f9p, -D-0.7c/0.57c -N
plot_blocks

gmt basemap -JQ88/34/5c -R$Region1 -BWSen -Ba10f2 -X+6.75c #-Bxya5f5g5
#gmt grd2cpt $EIGEN6C4_FIELD"_grd" $EGM_FIELD"_grd" -Chaxby -E30
gmt grdimage $EGM_FIELD"_grd" -I+d -E300 -Q #-Chaxby
#gmt colorbar -DJBC+h+w3.6c/0.2c+o0.c/0.9c -Bxaf -By+L"mGal"
echo "(b) EGM2008" | gmt text -F+cTL+f9p, -D-0.7c/0.57c -N
plot_blocks

gmt basemap -JQ88/34/5c -R$Region1 -BWSen -Ba10f2 -X-6.75c -Y-5.25c #-Bxya5f5g5
#gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5
#gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
gmt grdimage "residual_grd" -I+d -E300 -Q #-Chaxby
gmt colorbar -DJBC+h+w4.9c/0.2c+o3.21c/1.0c -Bxaf -By+L"mGal" --MAP_FRAME_PEN=0.6p
echo "(c) Differences" | gmt text -F+cTL+f9p, -D-0.7c/0.57c -N
#        plot_blocks

gmt basemap -JX4.69c/3.49c -R-15/15/0/1900 -BWSrt -Bxa10f5+l"Residuals (mGal)" -Bya+l"Counts" -X7.05 -Y0.11 --MAP_FRAME_PEN=1.p
gmt histogram "residual_gmt_format" -W1 -L1p -i2
#		gmt histogram "residual_gmt_format" -JX? -Bxa+l"mGal" -Bya200+l"Counts" -BWSrt -W1+b -L1p -i2
echo "(d) Histogram of the differences" | gmt text -F+cTL+f9p, -D-0.99c/0.64c -N
meanstd=$(python get_mean_std.py residual_gmt_format)
echo "${meanstd}" | gmt text -F+cTL+f9p, -D0.1c/-0.2c -N
#gmt basemap -JL88/34/30/43/? -R69/107/23/45 -BWSen -Bxya5f5g5
#gmt basemap -JM88/34/? -R69/107/23/45 -BwSEn -Baf5g5
#gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
#gmt grdimage  $ISOSTASY_EFFECT"_grd"  -E200 -t30 -Q -Cpanoply
#gmt colorbar -DJBC+h+w7c/0.2c+o0.c/1c -Bxa100f -By+L"mGal"
gmt end show
rm $EIGEN6C4_FIELD"_gmt_format"
rm $EIGEN6C4_FIELD"_grd"
rm $EGM_FIELD"_gmt_format"
rm $EGM_FIELD"_grd"
