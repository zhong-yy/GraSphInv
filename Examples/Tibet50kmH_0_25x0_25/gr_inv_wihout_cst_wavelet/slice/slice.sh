function check_and_delete(){
    if [ -f "$1" ];then
        echo "delete $1"
        rm $1
    fi
}

model=../gr_inv_no_constraint.txt
topo=../../tibet.nc
if [ ! -f "$topo" ]; then
    echo "Downloading topography data"
    gmt grdcut @earth_relief_15s.grd -R64/112/19/49 -G$topo
else
    echo "${topo} exists"
fi
#gmt grdcut @earth_relief_15s.grd -R69/107/23/45 -G$topo
#s0=$(seq 40 20 400)

#a0=($s0)
#len_a0=${#a0[*]}                          #a0的长度
#s1=$(seq ${a0[$len_a0 - 1]} -20 ${a0[0]}) #序列s0的倒序
#s1=$(seq 400 -40 40) #list

#a1=($s1) #array
#for i in $s1; do
#	rslice=$i.nc
#	r=$(echo "scale=5;6378137-$i*1000" | bc)
#	$extract_slice_tool $model $rslice r $r
#	echo "$i $r"
#done

start=('74.75/27.79' '80.70/25.48' '85.50/23.9' '91.86/21.78')
end=('78.20/40.0' '83.90/37.5' '88.5/36.2' '94.31/34.18')


R_interp=("" "")
img_Region=("" "")
ele_Region=("" "")

loop=$(seq 0 3)

gmt math -o0 -T40/400/10 --FORMAT_FLOAT_OUT=%20.3f T = zpoints.txt   

dep0=40.0
dep1=400.0
for j in $loop; do
	gmt project -C${start[$j]} -E${end[$j]} -G20 -Q -V >track$j.data
    
    check_and_delete cross_section$j.txt
    while read z; do
    	#gawk -v PREC=100 -v dep=${z} '{print $1, $2, $3, dep}' track$j.data  >> cross_section$j.txt
        gawk -v PREC=100 -v dep=${z} '{printf "%20.6f %20.6f %20.6f %20.6f\n", $1, $2, $3, dep}' track$j.data  >> cross_section$j.txt
    done < zpoints.txt
    python get_profile.py $model cross_section$j.txt track$j.profile

	gmt project -C${start[$j]} -E${end[$j]} -G10 -Q -V >track_e$j.data #for elevation
	R_interp[j]=$(gmt gmtinfo -i2 -C track$j.data | gawk '{print $1"/"$2}')"/$dep0/$dep1"
	img_Region[j]=$(gmt gmtinfo -i2 -C track$j.data | gawk '{print $1"/"$2}')"/0/400"
    
	gmt surface track$j.profile -R${R_interp[j]} -I50/20 -V -Gslice$j

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

gmt begin vslice_without_constraints jpg,eps E300
gmt set FONT 9p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_FRAME_PEN 1.2p,black
#gmt set MAP_FRAME_TYPE fancy+
#	gmt makecpt -Cwysiwyg -T-90/90/10 -D -Ic -H > mycpt.cpt
gmt makecpt -Cwysiwyg -T-90/90/10 -Ic -D -H >mycpt.cpt
#	gmt grd2cpt slice0 slice1 slice2 -Ic -D -H -Cwysiwyg -E30 >mycpt.cpt

#A
gmt basemap -Jx0.005/-0.005 -R${img_Region[0]} -Bxaf -Bya200f100g100+l"Depth (km)" -BWtsr
gmt grdimage slice0 -E300 -Q -Cmycpt.cpt
gmt plot moho0 -W0.5p,black,dashed -Gwhite -t30
gmt plot moho0 -Wthick,black,dashed

gmt plot LAB0 -Wthick,RED2,dashed
echo "(a)" | gmt text -F+cTL+f10p, -D-1.45c/0.5c -N
#gmt colorbar -DJBC+h+w5c/0.18c+o0.c/1.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt

Ra=($(gmt gmtinfo -I- -C elevation0))
ele0=${Ra[2]}
ele1=$(echo "${Ra[3]}+2.2" | bc)
gmt basemap -Jx0.005/0.06 -R${Ra[0]}/${Ra[1]}/$ele0/$ele1 -BWrb -Bya6f2 -Y2c
gmt plot elevation0 -W0.5p,black -Ggray

gmt basemap -Jx0.005/0.08 -R${Ra[0]}/${Ra[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 A
1350 7.6 9p 0 A\234
EOF
echo 451.82 4.8 -90 0.32 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 451.82 6.55 8p 0 MFT | gmt text -F+f+A
echo 800 8.3 8p 0 Tibetan plateau | gmt text -F+f+A
echo 1250 3.2 8p 0 Tarim | gmt text -F+f+A

#B
Rb=($(gmt gmtinfo -I- -C elevation1))
gmt basemap -Jx0.005/-0.005 -R${img_Region[1]} -Bxaf -Bya200f100g100+l"Depth (km)" -BWtsr -Y-5.25c
gmt grdimage slice1 -E300 -Q -Cmycpt.cpt
gmt plot moho1 -W0.5p,black,dashed -Gwhite -t30
gmt plot moho1 -Wthick,black,dashed

gmt plot LAB1 -Wthick,RED2,dashed
echo "(b)" | gmt text -F+cTL+f10p, -D-1.45c/0.5c -N
#gmt colorbar -DJBC+h+w5c/0.18c+o0.c/2.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt

gmt basemap -Jx0.005/0.06 -R${Rb[0]}/${Rb[1]}/$ele0/$ele1 -BWrb -Bya6f2 -Y2c
gmt plot elevation1 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Rb[0]}/${Rb[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 B
1330 7.6 9p 0 B\234
EOF
echo 323.2 4.9 -90 0.33 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 323.2 6.4 8p 0 MFT | gmt text -F+f+A
echo 594.7 7.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 594.7 9. 8p 0 YZS | gmt text -F+f+A

echo 817.0 7.4 -90 0.24 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 817.0 9. 8p 0 BNS | gmt text -F+f+A

echo 1097 7.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 1097 9. 8p 0 JRS | gmt text -F+f+A

#C
Rc=($(gmt gmtinfo -I- -C elevation2))
gmt basemap -Jx0.005/-0.005 -R${img_Region[2]} -Bxaf -Bya200f100g100+l"Depth (km)" -BWtsr -Y-5.25c
gmt grdimage slice2 -E300 -Q -Cmycpt.cpt
gmt plot moho2 -W0.5p,black,dashed -Gwhite -t30
gmt plot moho2 -Wthick,black,dashed

gmt plot LAB2 -Wthick,RED2,dashed
echo "(c)" | gmt text -F+cTL+f10p, -D-1.45c/0.5c -N
#gmt colorbar -DJBC+h+w5c/0.18c+o0.c/2.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt

gmt basemap -Jx0.005/0.06 -R${Rc[0]}/${Rc[1]}/$ele0/$ele1 -BWrb -Bya6f2 -Y2c
gmt plot elevation2 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Rc[0]}/${Rc[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 C
1350 7.6 9p 0 C\234
EOF
echo 372.70 5.0 -90 0.27 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 372.70 6.55 8p 0 MFT | gmt text -F+f+A
echo 606.63 7.5 -90 0.25 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 606.63 9.1 8p 0 YZS | gmt text -F+f+A
echo 931.74 7.5 -90 0.25 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 931.74 9.1 8p 0 BNS | gmt text -F+f+A

#D
gmt basemap -Jx0.005/-0.005 -R${img_Region[3]} -Bxaf+l"Distance (km)" -Bya200f100g100+l"Depth (km)" -BWtSr -Y-5.25c
gmt grdimage slice3 -E300 -Q -Cmycpt.cpt
gmt plot moho3 -W0.5p,black,dashed -Gwhite -t30
gmt plot moho3 -Wthick,black,dashed

gmt plot LAB3 -Wthick,RED2,dashed
echo "(d)" | gmt text -F+cTL+f10p, -D-1.45c/0.5c -N
#	gmt colorbar -DJBC+w5.5c/0.25c+o-3.2c/0.9c+e -Bxa20f10 -By+L"kg/m@+3@+" -C$mycpt
gmt colorbar -DJBC+w5c/0.2c+o0.c/1.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt

Rd=($(gmt gmtinfo -I- -C elevation3))
gmt basemap -Jx0.005/0.06 -R${Rd[0]}/${Rd[1]}/$ele0/$ele1 -BWrb -Bya6f2 -Y2c
gmt plot elevation3 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Rd[0]}/${Rd[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 D
1350 7.6 9p 0 D\234
EOF
echo 596.20 5.0 -90 0.35 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 596.20 6.55 8p 0 MFT | gmt text -F+f+A
echo 815.85 7.6 -90 0.25 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 815.85 9.2 8p 0 YZS | gmt text -F+f+A

echo 1129.64 7.6 -90 0.25 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 1129.64 9.2 8p 0 BNS | gmt text -F+f+A

#A2
gmt basemap -Jx0.005/-0.005 -R${img_Region[0]} -Bxaf -Bya200f100+l"Depth (km)" -Bltse -X7.2c -Y7.75c
Ra=($(gmt gmtinfo -I- -C elevation0))
ele0=${Ra[2]}
ele1=$(echo "${Ra[3]}+2.2" | bc)
gmt basemap -Jx0.005/0.06 -R${Ra[0]}/${Ra[1]}/$ele0/$ele1 -Blrb -Bya6f2 -Y2c
gmt plot elevation0 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Ra[0]}/${Ra[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 A
1350 7.6 9p 0 A\234
EOF
echo 451.82 4.8 -90 0.32 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 451.82 6.55 8p 0 MFT | gmt text -F+f+A
echo 800 8.3 8p 0 Tibetan plateau | gmt text -F+f+A
echo 1250 3.2 8p 0 Tarim | gmt text -F+f+A

#B2
Rb=($(gmt gmtinfo -I- -C elevation1))
gmt basemap -Jx0.005/-0.005 -R${img_Region[1]} -Bxaf -Bya200f100+l"Depth (km)" -Bltse -Y-5.25c
gmt basemap -Jx0.005/0.06 -R${Rb[0]}/${Rb[1]}/$ele0/$ele1 -Blrb -Bya6f2 -Y2c
gmt plot elevation1 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Rb[0]}/${Rb[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 B
1330 7.6 9p 0 B\234
EOF
echo 323.2 4.9 -90 0.33 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 323.2 6.4 8p 0 MFT | gmt text -F+f+A
echo 594.7 7.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 594.7 9. 8p 0 YZS | gmt text -F+f+A

echo 817.0 7.4 -90 0.24 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 817.0 9. 8p 0 BNS | gmt text -F+f+A

echo 1097 7.4 -90 0.22 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 1097 9. 8p 0 JRS | gmt text -F+f+A
#C2
Rc=($(gmt gmtinfo -I- -C elevation2))
gmt basemap -Jx0.005/-0.005 -R${img_Region[2]} -Bxaf -Bya200f100+l"Depth (km)" -Bltse -Y-5.25c
gmt basemap -Jx0.005/0.06 -R${Rc[0]}/${Rc[1]}/$ele0/$ele1 -Blrb -Bya6f2 -Y2c
gmt plot elevation2 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Rc[0]}/${Rc[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 C
1350 7.6 9p 0 C\234
EOF
echo 372.70 5.0 -90 0.27 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 372.70 6.55 8p 0 MFT | gmt text -F+f+A
echo 606.63 7.5 -90 0.25 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 606.63 9.1 8p 0 YZS | gmt text -F+f+A
echo 931.74 7.5 -90 0.25 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
echo 931.74 9.1 8p 0 BNS | gmt text -F+f+A
#D2
gmt basemap -Jx0.005/-0.005 -R${img_Region[3]} -Bxaf+l"Distance (km)" -Bya200f100+l"Depth (km)" -BltSe -Y-5.25c
Rd=($(gmt gmtinfo -I- -C elevation3))
gmt basemap -Jx0.005/0.06 -R${Rd[0]}/${Rd[1]}/$ele0/$ele1 -Blrb -Bya6f2 -Y2c
gmt plot elevation3 -W0.5p,black -Ggray
gmt basemap -Jx0.005/0.08 -R${Rd[0]}/${Rd[1]}/$ele0/$(echo "$ele1+2" | bc) -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 D
1350 7.6 9p 0 D\234
EOF
echo 596.20 5.0 -90 0.35 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 596.20 6.55 8p 0 MFT | gmt text -F+f+A
echo 815.85 7.6 -90 0.25 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 815.85 9.2 8p 0 YZS | gmt text -F+f+A

echo 1129.64 7.6 -90 0.25 | gmt plot -Sv0.15c+eA -W0.6p -Gblack
echo 1129.64 9.2 8p 0 BNS | gmt text -F+f+A

gmt end show


for j in $loop; do
rm moho$j.profile
rm elevation$j
rm slice$j
rm track_e$j.data
rm track_e$j.profile
rm moho$j
rm track$j.data
#rm track$j.profile
done
