model=../ref_model.nc
extract_slice_tool=Extract_slice
rslice1=./model_r_slice1
rslice2=./model_r_slice2
rslice3=./model_r_slice3
rslice4=./model_r_slice4
$extract_slice_tool $model $rslice1 r 6278137 ref
$extract_slice_tool $model $rslice2 r 6228137 ref
$extract_slice_tool $model $rslice3 r 6178137 ref
$extract_slice_tool $model $rslice4 r 6128137 ref

num_lon=185
num_lat=113
xyinc=${num_lon}+n/${num_lat}+n

Region1=65/111/20/48

gmt grdfilter $rslice1 -G${rslice1} -Fg160 -D4 -R${region} -I$xyinc -V
gmt grdfilter $rslice2 -G${rslice2} -Fg160 -D4 -R${region} -I$xyinc -V
gmt grdfilter $rslice3 -G${rslice3} -Fg160 -D4 -R${region} -I$xyinc -V
gmt grdfilter $rslice4 -G${rslice4} -Fg160 -D4 -R${region} -I$xyinc -V


function plot_blocks() {
        # ============绘制板块边界
        #gmt plot CN-plate-neighbor.dat -W1.0p,black -Sf0.5+t+l -G2/138/210
        #gmt plot Eurasia_India_PB.dat -W0.9p,black,-
        # ============绘制地块边界
        #gmt plot China_tectonic.dat -W0.3p,black,-
        # ============腾冲火山
        #echo 98 25 |gmt plot -St0.25 -Gwhite -W0.3p,black
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
function plot_profiles() {

	mark_start=('A' 'B' 'C' 'D')
	mark_end=("A'" "B'" "C'" "D'")

	mark_start=('A' 'B' 'C' 'D')
	mark_end=("A\234" "B\234" "C\234" "D\234")

	start=('74.75/27.79' '80.70/25.48' '85.50/23.9' '91.86/21.78')
	end=('78.20/40.0' '83.90/37.5' '88.5/36.2' '94.31/34.18')
	start_text=('73.65/27.5' '79.50/25.20' '84.10/23.80' '90.56/21.88')
	end_text=('76.45/39.50' '82.50/37.70' '87.30/36.5' '93.11/34.48')

	loop=$(seq 0 3)
	for j in $loop; do
		gmt project -C${start[$j]} -E${end[$j]} -G25 -Q -V >track$j.data
	done

	for j in $loop; do
		#gmt project -C${start[$j]} -E${end[$j]} -G25 | gmt plot -W1p,blue,
		gmt plot -W0.8p,black, track$j.data
	done

	rm annotations
	for i in $loop; do
		echo ${start_text[$i]}" 8.5p,black 0 "${mark_start[i]} | gawk '{gsub("/"," ");print}' >>annotations
		echo ${end_text[$i]}" 8.5p,black 0 "${mark_end[i]} | gawk '{gsub("/"," ");print}' >>annotations
	done
	gmt text -F+f+A annotations
}
gmt begin hslices_ini jpg,eps E300
gmt set FONT 9p
gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP ddd:mm:ssF
gmt set MAP_FRAME_PEN 1.5p,black

#origin of X/Y options is the upper left corner

intensity=+a-45+nt0.75+m0
contours=-45,-15,15,45+o+gGRAY+f7p
ct_lb_d=d5c
gmt makecpt -Cwysiwyg -T-60/60/0.1 -Z -Ic -H -D >mycpt.cpt

#		gmt basemap -JL88/34/34/45/4.8c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5
gmt basemap -JQ88/34/6c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5
gmt grdimage ${rslice1} -E300 -Q -nc  -Cmycpt.cpt #-I$intensity  #-I+a45+nt1+m0
#gmt grdcontour ${rslice1} -A$contours -G$ct_lb_d
plot_profiles
plot_blocks
#plot_ma
echo "(a) 100km depth" | gmt text -F+cTL+f9p, -D-0.9c/0.57c -N

#gmt colorbar -Cmycpt.cpt -D0.5/0.25+jCM+w7c/0.2c+h+e+n -B+l"Master CPT"

#-JL88/34/34/45/4.8c
gmt basemap -JQ88/34/6c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5 -X7.55c
#gmt grdimage @earth_relief_01m -Cetopo1 -t20 -I+d
gmt grdimage ${rslice2} -E300 -Q -nc  -Cmycpt.cpt #-I$intensity  #-I+a45+nt1+m0
#gmt grdcontour ${rslice2} -A$contours -G$ct_lb_d
plot_profiles
plot_blocks
#plot_ma
echo "(b) 150km depth" | gmt text -F+cTL+f9p, -D-0.9c/0.57c -N

gmt basemap -JQ88/34/6c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5 -X-7.55c -Y-5.75c
#gmt grdimage @earth_relief_01m -Cetopo1 -t20 -I+d
gmt grdimage ${rslice3} -E300 -Q -nb  -Cmycpt.cpt #-I$intensity  #-I+a45+nt1+m0
#gmt grdcontour ${rslice3} -A$contours -G$ct_lb_d
plot_profiles
plot_blocks
#plot_ma
echo "(c) 200km depth" | gmt text -F+cTL+f9p, -D-0.85c/0.57c -N

gmt basemap -JQ88/34/6c -R$Region1 -BWSen -Bxa10f2g10 -Bya5f1g5 -X7.55c
#gmt grdimage @earth_relief_01m -Cetopo1 -t20 -I+d
gmt grdimage ${rslice4} -E300 -Q -nc  -Cmycpt.cpt #-I$intensity  #-I+a45+nt1+m0
#gmt grdcontour ${rslice4} -A$contours -G$ct_lb_d
#		gmt colorbar -DJBC+h+w4.5c/0.2c+o0.c/0.8c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt
gmt colorbar -Dn-0.6/-0.25+h+w5.3c/0.2c+o0.c/-0.c+e -Bxa50f10 -By+L"kg/m@+3@+" -Cmycpt.cpt --MAP_FRAME_PEN=0.6p
plot_profiles
plot_blocks
#plot_ma
echo "(d) 250km depth" | gmt text -F+cTL+f9p, -D-0.9c/0.57c -N
gmt end show

rm $rslice1
rm $rslice2
rm $rslice3
rm $rslice4
for j in $loop; do
	rm track$j.data
done
rm annotations
