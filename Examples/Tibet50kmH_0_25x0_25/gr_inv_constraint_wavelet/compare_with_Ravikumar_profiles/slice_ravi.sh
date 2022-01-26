function check_and_delete(){
    if [ -f "$1" ];then
        echo "delete $1"
        rm $1
    fi
}
model=../gr_inv_pet.txt
#extract_slice_tool=Extract_slice
topo=../../tibet.nc
if [ ! -f "$topo" ]; then
    echo "Downloading topography data"
    gmt grdcut @earth_relief_15s.grd -R64/112/19/49 -G$topo
else
    echo "${topo} exists"
fi
dep0=40.0
dep1=400.0

start=('76.24/28.14' '84.86/26.42' '87.33/25.51')
end=('78.96/39.85' '84.84/39.33' '97.23/39.64')

mark_start=('A' 'B' 'C')
mark_end=("A'" "B'" "C'")

R_interp=("" "")
img_Region=("" "")
ele_Region=("" "")

loop=$(seq 0 2)

gmt math -o0 -T40/400/10 --FORMAT_FLOAT_OUT=%20.3f T = zpoints.txt   
for j in $loop; do
	gmt project -C${start[$j]} -E${end[$j]} -G20 -Q -V >track$j.data
    check_and_delete cross_section$j.txt
    while read z; do
    	#gawk -v PREC=100 -v dep=${z} '{print $1, $2, $3, dep}' track$j.data  >> cross_section$j.txt
        gawk -v PREC=100 -v dep=${z} '{printf "%20.6f %20.6f %20.6f %20.6f\n", $1, $2, $3, dep}' track$j.data  >> cross_section$j.txt
    done < zpoints.txt
    python get_profile.py $model cross_section$j.txt track$j.profile

	gmt project -C${start[$j]} -E${end[$j]} -G10 -Q -V >track_e$j.data #for elevation"
	R_interp[j]=$(gmt gmtinfo -i2 -C track$j.data | gawk '{print $1"/"$2}')"/$dep0/$dep1"
	img_Region[j]=$(gmt gmtinfo -i2 -C track$j.data | gawk '{print $1"/"$2}')"/0/400"
	gmt surface track$j.profile -R${R_interp[j]} -I25/20 -V -Gslice$j
	#echo "gmt surface track$j.profile -R${R_interp[j]} -I100/50 -V -Gslice$j"

	#extract elevation
	gmt grdtrack track_e$j.data -G$topo >track_e$j.profile

	#tac命令将文件逆序输出
	tac track_e$j.data | gawk '{print $1,$2,$3,0}' >>track_e$j.profile

	gawk '{printf"%15.6f %15.6f\n",$3,$4/1000.0}' track_e$j.profile >elevation$j

	#Extract depth to the moho along profiles, with respect to sea level
	gmt grdtrack track$j.data -Gdepthtomoho.nc >moho$j.profile
	tac moho$j.profile | gawk '{print $1,$2,$3,0}' >>moho$j.profile
	gawk '{printf"%15.6f %15.6f\n",$3,-1.0*$4}' moho$j.profile >moho$j
done

contours=10,20,30+o+gGRAY+f7p
contours2=-30,-20,-10+o+gGRAY+f7p
ct_lb_d=d1.5c
gmt begin slices_Ravikumar jpg,eps E300
gmt set FONT 9p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_FRAME_PEN 1.2p,black
#gmt set MAP_FRAME_TYPE fancy+
#	gmt makecpt -Cwysiwyg -T-90/90/10 -D -Ic -H > mycpt.cpt
gmt makecpt -Cwysiwyg -T-80/80/0.1 -Z -Ic -D -H >mycpt.cpt
#	gmt grd2cpt slice0 slice1 slice2 -Ic -D -H -Cwysiwyg -E30 >mycpt.cpt

#A
gmt basemap -Jx0.005/-0.005 -R${img_Region[0]} -Bxaf -Bya200f100g100+l"Depth (km)" -BWtSr
gmt grdimage slice0 -E300 -Q -Cmycpt.cpt
gmt plot moho0 -W0.5p,black,dashed -Gwhite -t30
gmt plot moho0 -Wthick,black,dashed

gmt plot LAB0 -Wthick,RED2,dashed
gmt plot ./Ravi_lab/LAB_profileA.txt -Wthick,PURPLE3,dashed

echo "(a)" | gmt text -F+cTL+f10p, -D-1.45c/0.5c -N
#gmt colorbar -DJBC+h+w5c/0.18c+o0.c/1.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt



Ra=($(gmt gmtinfo -I- -C elevation0))
ele0=${Ra[2]}
ele1=$(echo "${Ra[3]}+4.2" | bc)
gmt basemap -Jx0.005/0.06 -R${Ra[0]}/${Ra[1]}/$ele0/$ele1 -BWrb -Bya6f2 -Y2c
gmt plot elevation0 -W0.5p,black -Ggray


gmt basemap -Jx0.005/0.08 -R${Ra[0]}/${Ra[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
35 7.6 9p 1 E
1290 7.6 9p 1 E\234
EOF
echo 290.6 4.8 -90 0.32 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 290.6 6.4 8p 0 MFT | gmt text -F+f+A
#echo 648.53 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 648.53 8. 8p 90 YZS | gmt text -F+f+A
#echo 713.19 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 713.19 8. 8p 90 BNS | gmt text -F+f+A
echo 871.08 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 871.08 8. 8p 0 JRS | gmt text -F+f+A
#echo 700 7.8 8p 0 Tibetan plateau | gmt text -F+f+A
echo 1150 3.2 8p 0 Tarim | gmt text -F+f+A

#B
Rb=($(gmt gmtinfo -I- -C elevation1))
gmt basemap -Jx0.005/-0.005 -R${img_Region[1]} -Bxaf -Bya200f100g100+l"Depth (km)" -BWtSr -Y-5.5c
gmt grdimage slice1 -E300 -Q -Cmycpt.cpt
gmt plot moho1 -W0.5p,black,dashed -Gwhite -t30
gmt plot moho1 -Wthick,black,dashed


gmt plot LAB1 -Wthick,RED2,dashed
gmt plot ./Ravi_lab/LAB_profileB.txt -Wthick,PURPLE3,dashed
echo "(b)" | gmt text -F+cTL+f10p, -D-1.45c/0.5c -N
#gmt colorbar -DJBC+h+w5c/0.18c+o0.c/2.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt

gmt colorbar -DJMR+w4c/0.22c+o0.75c/1.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt

gmt basemap -Jx0.005/0.06 -R${Rb[0]}/${Rb[1]}/$ele0/$ele1 -BWrb -Bya6f2 -Y2c
gmt plot elevation1 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Rb[0]}/${Rb[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
35 7.6 9p 0 F
1380 7.6 9p 0 F\234
EOF
echo 98.35 2.8 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 98.35 4.5 8p 0 MFT | gmt text -F+f+A

echo 338.84 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 338.84 8. 8p 0 YZS | gmt text -F+f+A

echo 641.27 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 641.27 8. 8p 0 BNS | gmt text -F+f+A

echo 958.25 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 958.25 8. 8p 0 JRS | gmt text -F+f+A

#echo 1186 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
#echo 1186 10. 6.5p 0 Altyn Tagh | gmt text -F+f+A
#echo 1186 7.6 6.5p 0 Fault | gmt text -F+f+A
#C
Rc=($(gmt gmtinfo -I- -C elevation2))
gmt basemap -Jx0.005/-0.005 -R${img_Region[2]} -Bxaf+l"Distance (km)" -Bya200f100g100+l"Depth (km)" -BWtSr -Y-5.5c
gmt grdimage slice2 -E300 -Q -Cmycpt.cpt
gmt plot moho2 -W0.5p,black,dashed -Gwhite -t30
gmt plot moho2 -Wthick,black,dashed
gmt plot ./Ravi_lab/LAB_profileC.txt -Wthick,PURPLE3,dashed

gmt plot LAB2 -Wthick,RED2,dashed
echo "(c)" | gmt text -F+cTL+f10p, -D-1.45c/0.5c -N

#gmt colorbar -DJBC+h+w5c/0.18c+o0.c/2.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt

gmt legend  -DjBC+w8.8c+o0.c/-1.5c --FONT_ANNOT_PRIMARY=8p <<EOF
N 3.8c 5.c
S 0.15c - 0.5 - thin,RED2,- 0.65c LAB from LITHO1.0
S 0.15c - 0.5 - thin,PURPLE3,- 0.65c LAB from Ravikumar et al., 2020
EOF

gmt basemap -Jx0.005/0.06 -R${Rc[0]}/${Rc[1]}/$ele0/$ele1 -BWrb -Bya6f2 -Y2c
gmt plot elevation2 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Rc[0]}/${Rc[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
35 7.6 9p 0 G
1780 7.6 9p 0 G\234
EOF
echo 152.37 2.8 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 152.37 4.5 8p 0 MFT | gmt text -F+f+A

echo 467.31 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 467.31 8. 8p 0 YZS | gmt text -F+f+A

echo 846.69 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 846.69 8. 8p 0 BNS | gmt text -F+f+A

echo 1134.85 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 1134.85 8. 8p 0 JRS | gmt text -F+f+A

echo 1298 6.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 1298 10. 6.5p 0 Kunlun | gmt text -F+f+A
echo 1298 7.6 6.5p 0 Fault | gmt text -F+f+A
gmt end show

for j in $loop; do
rm moho$j.profile
rm elevation$j
rm slice$j
rm track_e$j.data
rm track_e$j.profile
rm moho$j
rm track$j.data
rm track$j.profile
done
