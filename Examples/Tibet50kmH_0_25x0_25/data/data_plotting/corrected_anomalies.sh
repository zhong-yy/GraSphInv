TOTAL=../Eigen6C4/Tibet_total_g_r
NORMAL=../Eigen6C4/Normal_g_r
FREEAIR=../Eigen6C4/Free_air_g_r
BOUGUER=../Eigen6C4/Bouguer_g_r
RM_SEDIM_CRYSTL=../Crust_Correction/rm_sedim_cryst_g_r
MANTLE=../Crust_Correction/mantle_g_r

echo $OBSERVED_FIELD
echo $PREDICTED_FIELD

python swap_2cols_of_g_data.py $FREEAIR "_gmt_format"
python swap_2cols_of_g_data.py $TOTAL "_gmt_format"
python swap_2cols_of_g_data.py $NORMAL "_gmt_format"
python swap_2cols_of_g_data.py $BOUGUER "_gmt_format"
python swap_2cols_of_g_data.py $RM_SEDIM_CRYSTL "_gmt_format"
python swap_2cols_of_g_data.py $MANTLE "_gmt_format"

Region1=65/111/20/48
xinc_yinc=0.25/0.25
gmt xyz2grd $FREEAIR"_gmt_format" -I$xinc_yinc -R$Region1 -G$FREEAIR"_grd"
gmt xyz2grd $TOTAL"_gmt_format" -I$xinc_yinc -R$Region1 -G$TOTAL"_grd"
gmt xyz2grd $NORMAL"_gmt_format" -I$xinc_yinc -R$Region1 -G$NORMAL"_grd"
gmt xyz2grd $BOUGUER"_gmt_format" -I$xinc_yinc -R$Region1 -G$BOUGUER"_grd"
gmt xyz2grd $RM_SEDIM_CRYSTL"_gmt_format" -I$xinc_yinc -R$Region1 -G$RM_SEDIM_CRYSTL"_grd"
gmt xyz2grd $MANTLE"_gmt_format" -I$xinc_yinc -R$Region1 -G$MANTLE"_grd"

#lon_inv=(67.5 109.5)
#lat_inv=(21.5 47.5)

#inversion_region='inversion_region'

#python make_frame.py ${lon_inv[0]} ${lon_inv[1]} ${lat_inv[0]} ${lat_inv[1]} $inversion_region

function plot_blocks() {
# ============Main frontal thrust
#gmt plot Eurasia_India_PB.dat -W1.0p,red,-
gmt plot MFT.gmt  -Wthinner,black -Sf1.c/0.1c+l+t+o0.2c -Gblack

# ============suture
gmt plot sutures.gmt -Wthinner,black,-

# ============thrust belt
#gmt plot western_kunlun.gmt -Wthinner,black, -Sf2c/0.3c+l+s+o1
#gmt plot qilian_shan.gmt -Wthinner,black, #-Sf2c/0.3c+l+s+o1c

#=============Fault
gmt plot kunlun_fault.gmt -Wthinner,black, #-Sf2c/0.3c+l+s+o0.8c
gmt plot Altyn_Tagh_Fault.gmt -Wthinner,black, -Sf2c/0.3c+l+s+o1c
gmt plot Karakax.gmt -Wthinner,black, #-Sf2c/0.3c+l+s+o1c
gmt plot Karakoram_Fault.gmt -Wthinner,black, #-Sf2c/0.3c+l+s+o1c
gmt plot Haiyuan_Fault.gmt -Wthinner,black, #-Sf2c/0.3c+l+s+o1c
gmt plot Xianshuihe_Fault.gmt -Wthinner,black, #-Sf2c/0.3c+l+s+o1c
gmt plot Sagaing_Fault.gmt -Wthinner,black, #-Sf2c/0.3c+l+s+o1c
gmt plot Red_river_fault.gmt -Wthinner,black, #-Sf2c/0.3c+l+s+o1c
        # ============绘制板块边界
        #gmt plot CN-plate-neighbor.dat -W1.0p,black -Sf0.5+t+l -G2/138/210
        #gmt plot Eurasia_India_PB.dat -W0.9p,black,-
        # ============绘制地块边界
        #gmt plot China_tectonic.dat -W0.3p,black,-
        # ============腾冲火山
        #echo 98 25 |gmt plot -St0.25 -Gwhite -W0.3p,black
}
CPT=haxby
gmt begin gr_ano jpg,eps E300
gmt set FONT_TAG 9p
gmt set FONT_TITLE 9p
gmt set MAP_TITLE_OFFSET 9p
gmt set FONT 9p
gmt set MAP_FRAME_TYPE plain
#	gmt set MAP_FRAME_WIDTH 2.5p
gmt set MAP_FRAME_PEN 1.5p
gmt set FORMAT_GEO_MAP ddd:mm:ssF
#gmt basemap -JL88/34/30/43/4.5c -R$Region1 -BWSen -Bafg #-Bxya5f5g5
gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Bxya10f2
gmt grd2cpt $TOTAL"_grd" -C${CPT} -Ic -E -Z
gmt grdimage $TOTAL"_grd" -E300 -Q -I+d #-C${CPT}
#gmt grdcontour $TOTAL"_grd" -C100 -A300
plot_blocks
gmt colorbar -DJBC+h+w3.6/0.18c+o0.c/0.8c -Bxa500f -By+L"mGal" --MAP_FRAME_PEN=0.6p
echo "(a)" | gmt text -F+cTL+f11p, -D-0.8c/0.15c -N

#gmt basemap -JL88/34/30/43/4.5c -R$Region1 -BWSen -Bafg -X+6.4c #-Bxya5f5g5
gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Ba10f2g -X+6.5c
gmt grd2cpt $NORMAL"_grd" -C${CPT} -Ic -E -Z
gmt grdimage $NORMAL"_grd" -E300 -Q  -I+d #-C${CPT}
gmt colorbar -DJBC+h+w3.6/0.18c+o0.c/0.8c -Bxa500f -By+L"mGal" --MAP_FRAME_PEN=0.6p
echo "(b)" | gmt text -F+cTL+f11p, -D-0.8c/0.15c -N
plot_blocks

gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Ba10f2g -X-6.5c -Y-5.5c #-Bxya5f5g5
#gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5
#gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
gmt grd2cpt $FREEAIR"_grd" -C${CPT} -Ic -E -Z
gmt grdimage $FREEAIR"_grd" -E300 -Q  -I+d #-C${CPT}
plot_blocks
gmt colorbar -DJBC+h+w3.6/0.18c+o0.c/0.8c -Bxaf -By+L"mGal" --MAP_FRAME_PEN=0.6p
echo "(c)" | gmt text -F+cTL+f11p, -D-0.8c/0.15c -N

gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Ba10f2g -X+6.5c #-Bxya5f5g5
#gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5
#gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
gmt grd2cpt $BOUGUER"_grd" -C${CPT} -Ic -E -Z
gmt grdimage $BOUGUER"_grd" -E300 -Q  -I+d #-C${CPT}
plot_blocks
gmt colorbar -DJBC+h+w3.6/0.18c+o0.c/0.8c -Bxaf -By+L"mGal" --MAP_FRAME_PEN=0.6p
echo "(d)" | gmt text -F+cTL+f11p, -D-0.8c/0.15c -N

gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Ba10f2g -X-6.5c -Y-5.5c #-Bxya5f5g5
#gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5
#gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
gmt grd2cpt $RM_SEDIM_CRYSTL"_grd" -C${CPT} -Ic -E -Z
gmt grdimage $RM_SEDIM_CRYSTL"_grd" -E300 -Q  -I+d #-C${CPT}
plot_blocks
gmt colorbar -DJBC+h+w3.6/0.18c+o0.c/0.8c -Bxaf -By+L"mGal" --MAP_FRAME_PEN=0.6p
echo "(e)" | gmt text -F+cTL+f11p, -D-0.8c/0.15c -N

gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Ba10f2g -X+6.5c #-Bxya5f5g5
#gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5
#gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
gmt grd2cpt $MANTLE"_grd" -C${CPT} -Ic -E -Z
gmt grdimage $MANTLE"_grd" -E300 -Q  -I+d #-C${CPT}
plot_blocks
gmt colorbar -DJBC+h+w3.6/0.18c+o0.c/0.8c -Bxaf -By+L"mGal" --MAP_FRAME_PEN=0.6p
echo "(f)" | gmt text -F+cTL+f11p, -D-0.8c/0.15c -N
gmt end show
rm $FREEAIR"_gmt_format"
rm $FREEAIR"_grd"
rm $BOUGUER"_gmt_format"
rm $BOUGUER"_grd"
rm $RM_SEDIM_CRYSTL"_gmt_format"
rm $RM_SEDIM_CRYSTL"_grd"
rm $MANTLE"_gmt_format"
rm $MANTLE"_grd"
rm $TOTAL"_gmt_format"
rm $TOTAL"_grd"
rm $NORMAL"_gmt_format"
rm $NORMAL"_grd"
