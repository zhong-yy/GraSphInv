model=../gr_inv_no_constraint.txt
extract_slice_tool=Extract_slice
rslice1=./model_r_slice1
rslice2=./model_r_slice2
rslice3=./model_r_slice3
rslice4=./model_r_slice4

start_lon=65
stop_lon=111
num_lon=185
start_lat=20
stop_lat=48
num_lat=113
python get_h_slice.py $model 6278.137 ${start_lon} ${stop_lon} ${num_lon} ${start_lat} ${stop_lat} ${num_lat} ${rslice1}.txt
python get_h_slice.py $model 6228.137 ${start_lon} ${stop_lon} ${num_lon} ${start_lat} ${stop_lat} ${num_lat} ${rslice2}.txt
python get_h_slice.py $model 6178.137 ${start_lon} ${stop_lon} ${num_lon} ${start_lat} ${stop_lat} ${num_lat} ${rslice3}.txt
python get_h_slice.py $model 6128.137 ${start_lon} ${stop_lon} ${num_lon} ${start_lat} ${stop_lat} ${num_lat} ${rslice4}.txt

xyinc=${num_lon}+n/${num_lat}+n
#$extract_slice_tool $model $rslice1 r 6278137
#$extract_slice_tool $model $rslice2 r 6228137
#$extract_slice_tool $model $rslice3 r 6178137
#$extract_slice_tool $model $rslice4 r 6128137

Region1=65/111/20/48
#gmt xyz2grd ${rslice1}.txt -I$xyinc -R$Region1 -G$rslice1
#gmt xyz2grd ${rslice2}.txt -I$xyinc -R$Region1 -G$rslice2
#gmt xyz2grd ${rslice3}.txt -I$xyinc -R$Region1 -G$rslice3
#gmt xyz2grd ${rslice4}.txt -I$xyinc -R$Region1 -G$rslice4

gmt surface ${rslice1}.txt -I$xyinc -R$Region1 -G$rslice1
gmt surface ${rslice2}.txt -I$xyinc -R$Region1 -G$rslice2
gmt surface ${rslice3}.txt -I$xyinc -R$Region1 -G$rslice3
gmt surface ${rslice4}.txt -I$xyinc -R$Region1 -G$rslice4

gmt grdfilter $rslice1 -G${rslice1} -Fg160 -D4 -R${region} -I$xyinc -V
gmt grdfilter $rslice2 -G${rslice2} -Fg160 -D4 -R${region} -I$xyinc -V
gmt grdfilter $rslice3 -G${rslice3} -Fg160 -D4 -R${region} -I$xyinc -V
gmt grdfilter $rslice4 -G${rslice4} -Fg160 -D4 -R${region} -I$xyinc -V

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

function plot_ma(){
gmt plot -Sc0.08 ./magmatism/Potassic.txt -Wthinnest,BLACK -GMAGENTA
gmt plot -Sc0.08 ./magmatism/Ultrapotasic_adakitc.txt -Wthinnest,BLACK -GORANGE2
}

function plot_a() {
	cat >annotations <<EOF
85.1 31.2 8.5p, -8 LB
85 33.5 8.5p, -8 QB
85 28.5 8.5p, -19 Himalayas
99.75 33.8 8.5p, -18.5 SGB
103.5 34.7 8.5p -13 Qinling Orogen
100.5 37.3 8.5p -15 Qilian Orogen
94.5 37.3 8.5p, -10 Qaidam Basin
85 39.5 8.5p, 0 Tarim Basin
81 42.8 8.5p, 10 Tien Shan
106 30.5 8.5p 10 Sichuan
105.5 29.5 8.5p 10 Basin
80.0 23.0 8.5p 0 Indian Plate
91.5 29.6 8.5p 0 YZS
92.0 32.3 8.5p 0 BNS
91.9 34.6 8.5p -10 JRS
84 27 8.5p -19 MFT
EOF
	gmt text -F+f+A annotations
}
gmt begin hslices_no_cst jpg,eps E300
gmt set FONT 9p
gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP ddd:mm:ssF
gmt set MAP_FRAME_PEN 1.5p,black

#origin of X/Y options is the upper left corner

intensity=+a-45+nt0.75+m0
contours=-45,-15,15,45+o+gGRAY+f7p
ct_lb_d=d5c
gmt makecpt -Cwysiwyg -T-80/80/0.1 -Z -Ic -H -D >mycpt.cpt

#		gmt basemap -JL88/34/34/45/4.8c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5
gmt basemap -JQ88/34/6c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5
gmt grdimage ${rslice1} -E300 -Q -nc -I$intensity -Cmycpt.cpt  #-I+a45+nt1+m0
#gmt grdcontour ${rslice1} -A$contours -G$ct_lb_d
#plot_profiles
plot_blocks
plot_ma
echo "(a) 100km depth" | gmt text -F+cTL+f9p, -D-0.9c/0.57c -N

#gmt colorbar -Cmycpt.cpt -D0.5/0.25+jCM+w7c/0.2c+h+e+n -B+l"Master CPT"

#-JL88/34/34/45/4.8c
gmt basemap -JQ88/34/6c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5 -X7.55c
#gmt grdimage @earth_relief_01m -Cetopo1 -t20 -I+d
gmt grdimage ${rslice2} -E300 -Q -nc -I$intensity -Cmycpt.cpt  #-I+a45+nt1+m0
#gmt grdcontour ${rslice2} -A$contours -G$ct_lb_d
#plot_profiles
plot_blocks
plot_ma
echo "(b) 150km depth" | gmt text -F+cTL+f9p, -D-0.9c/0.57c -N

gmt basemap -JQ88/34/6c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5 -X-7.55c -Y-5.75c
#gmt grdimage @earth_relief_01m -Cetopo1 -t20 -I+d
gmt grdimage ${rslice3} -E300 -Q -nb -I$intensity -Cmycpt.cpt  #-I+a45+nt1+m0
#gmt grdcontour ${rslice3} -A$contours -G$ct_lb_d
#plot_profiles
plot_blocks
plot_ma
echo "(c) 200km depth" | gmt text -F+cTL+f9p, -D-0.85c/0.57c -N

gmt basemap -JQ88/34/6c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5 -X7.55c
#gmt grdimage @earth_relief_01m -Cetopo1 -t20 -I+d
gmt grdimage ${rslice4} -E300 -Q -nc -I$intensity -Cmycpt.cpt  #-I+a45+nt1+m0
#gmt grdcontour ${rslice4} -A$contours -G$ct_lb_d
#		gmt colorbar -DJBC+h+w4.5c/0.2c+o0.c/0.8c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt
gmt colorbar -Dn-0.6/-0.25+h+w5.3c/0.2c+o0.c/-0.c+e -Bxa50f10 -By+L"kg/m@+3@+" -Cmycpt.cpt --MAP_FRAME_PEN=0.6p
#plot_profiles
plot_blocks
plot_ma
echo "(d) 250km depth" | gmt text -F+cTL+f9p, -D-0.9c/0.57c -N
gmt end show

rm $rslice1
rm $rslice2
rm $rslice3
rm $rslice4
rm ${rslice1}.txt
rm ${rslice2}.txt
rm ${rslice3}.txt
rm ${rslice4}.txt
