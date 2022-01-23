TOTAL=../Eigen6C4/Tibet_total_g_r
NORMAL=../Eigen6C4/Normal_g_r
TOPO=../Eigen6C4/Topo_effect_g_r
SEDIMENT=../Crust_Correction/sediments_g_r
CRYSTALLINE_CRUST=../Crust_Correction/crystalline_crust_g_r
MOHO=../Crust_Correction/moho_g_r

echo $OBSERVED_FIELD
echo $PREDICTED_FIELD

python swap_2cols_of_g_data.py $TOPO "_gmt_format"
python swap_2cols_of_g_data.py $TOTAL "_gmt_format"
python swap_2cols_of_g_data.py $NORMAL "_gmt_format"
python swap_2cols_of_g_data.py $SEDIMENT "_gmt_format"
python swap_2cols_of_g_data.py $CRYSTALLINE_CRUST "_gmt_format"
python swap_2cols_of_g_data.py $MOHO "_gmt_format"

Region1=65/111/20/48
xinc_yinc=0.25/0.25
gmt xyz2grd $TOPO"_gmt_format" -I$xinc_yinc -R$Region1 -G$TOPO"_grd"
gmt xyz2grd $TOTAL"_gmt_format" -I$xinc_yinc -R$Region1 -G$TOTAL"_grd"
gmt xyz2grd $NORMAL"_gmt_format" -I$xinc_yinc -R$Region1 -G$NORMAL"_grd"
gmt xyz2grd $SEDIMENT"_gmt_format" -I$xinc_yinc -R$Region1 -G$SEDIMENT"_grd"
gmt xyz2grd $CRYSTALLINE_CRUST"_gmt_format" -I$xinc_yinc -R$Region1 -G$CRYSTALLINE_CRUST"_grd"
gmt xyz2grd $MOHO"_gmt_format" -I$xinc_yinc -R$Region1 -G$MOHO"_grd"

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
}

gmt begin crustal_effects jpg,eps E300
	gmt set FONT_TAG 9p
	gmt set FONT_TITLE 9p
	gmt set MAP_TITLE_OFFSET 9p 
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 0.75p
	gmt set FORMAT_GEO_MAP ddd:mm:ssF
        gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Bafg #-Bxya5f5g5
        gmt grd2cpt $TOPO"_grd" -Chaxby -Ic  -Z -E300
        gmt grdimage  $TOPO"_grd"  -E300 -Q -I+d #-Chaxby
        gmt colorbar -DJBC+h+w3.6c/0.18c+o0.c/0.8c -Bxaf -By+L"mGal"
        echo "(a)" | gmt text -F+cTL+f9p, -D-0.65c/0.25c -N
        plot_blocks
        
        gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Bafg -X+6.4c #-Bxya5f5g5
        gmt grd2cpt $SEDIMENT"_grd" -Chaxby -Ic  -Z -E300
        gmt grdimage  $SEDIMENT"_grd"  -E300 -Q -I+d #-Chaxby
        gmt colorbar -DJBC+h+w3.6c/0.18c+o0.c/0.8c -Bxaf -By+L"mGal"
        echo "(b)" | gmt text -F+cTL+f9p, -D-0.65c/0.25c -N
        plot_blocks
        
        gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Bafg -X-6.4c -Y-5.4c #-Bxya5f5g5
        #gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5	        
        #gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
        gmt grd2cpt $CRYSTALLINE_CRUST"_grd" -Chaxby -Ic  -Z -E300
        gmt grdimage  $CRYSTALLINE_CRUST"_grd" -I+d  -E300 -Q #-Chaxby
        gmt colorbar -DJBC+h+w3.6c/0.18c+o0.c/0.8c -Bxaf -By+L"mGal"
        echo "(c)" | gmt text -F+cTL+f9p, -D-0.65c/0.25c -N
        plot_blocks
        
        gmt basemap -JQ88/34/5.0c -R$Region1 -BWSen -Bafg -X+6.4c #-Bxya5f5g5
        #gmt basemap -JM88/34/? -R69/107/23/45 -BWSen -Baf5g5	        
        #gmt grdimage @earth_relief_01m -R63/113/18/49 -Cetopo1 -t10 -I+
        gmt grd2cpt $MOHO"_grd" -Chaxby -Ic -Z -E300
        gmt grdimage  $MOHO"_grd"  -Q -I+d  -E300 #-Chaxby
        gmt colorbar -DJBC+h+w3.6c/0.18c+o0.c/0.8c -Bxaf -By+L"mGal"
        echo "(d)" | gmt text -F+cTL+f9p, -D-0.65c/0.25c -N
        plot_blocks
gmt end show
rm $TOPO"_gmt_format"
rm $TOPO"_grd"
rm $SEDIMENT"_gmt_format"
rm $SEDIMENT"_grd"
rm $CRYSTALLINE_CRUST"_gmt_format"
rm $CRYSTALLINE_CRUST"_grd"
rm $MOHO"_gmt_format"
rm $MOHO"_grd"
rm $TOTAL"_gmt_format"
rm $TOTAL"_grd"
rm $NORMAL"_gmt_format"
rm $NORMAL"_grd"
