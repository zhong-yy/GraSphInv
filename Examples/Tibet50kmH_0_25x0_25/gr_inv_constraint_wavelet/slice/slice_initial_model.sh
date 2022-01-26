model=../ref_model.nc
extract_slice_tool=Extract_slice
topo=../../tibet.nc
if [ ! -f "$topo" ]; then
    echo "Downloading topography data"
    gmt grdcut @earth_relief_15s.grd -R64/112/19/49 -G$topo
else
    echo "${topo} exists"
fi
s0=$(seq 40 40 400)

a0=($s0)
len_a0=${#a0[*]} #a0的长度
s1=$(seq ${a0[$len_a0-1]} -40 ${a0[0]}) #序列s0的倒序
#s1=$(seq 400 -40 40) #list

a1=($s1) #array
for i in $s1
do
rslice=$i.nc
r=`echo "scale=5;6378137-$i*1000"|bc`
$extract_slice_tool $model $rslice r $r ref
echo "$i $r"
done

dep0=40.0
dep1=400.0
r0=`echo "scale=3;(6378137.0-$dep1*1000.0)/1000.0"|bc`
r1=`echo "scale=3;(6378137.0-$dep0*1000.0)/1000.0"|bc`

#start=('86.00/23.50' '100.00/23.50' '85.0/34.5')    
#end=('90.00/35.0' '89.50/31.0' '76.5/44.5')

#start=('74.00/26' '80.0/24.0' '76.0/44')    
#end=('78.00/40.0' '84.0/38.0' '92.0/36.5')

start=('74.75/27.79' '80.70/25.48' '85.50/23.9' '91.86/21.78')
end=('78.20/40.0' '83.90/37.5' '88.5/36.2' '94.31/34.18')


mark_start=('A' 'B' 'C' 'D')
mark_end=("A'" "B'" "C'" "D'")


R_interp_polar=("" "")
img_Region_polar=("" "")
ele_Region_polar=("" "")

R_interp_cartesian=("" "")
img_Region_cartesian=("" "")
ele_Region_cartesian=("" "")

loop=$(seq 0 3)

#gmt xyz2grd ../../depthtomoho.xyz -Rg -I1/1 -Gout.grd -r -V

for j in $loop
do
    gmt project -C${start[$j]} -E${end[$j]} -G25 -Q -V > track$j.data
    gmt project -C${start[$j]} -E${end[$j]} -G10 -Q -V > track_e$j.data  #for elevation
    R_interp_cartesian[j]=`gmt gmtinfo -i2 -C track$j.data | gawk '{print $1"/"$2}'`"/$dep0/$dep1"
    img_Region_cartesian[j]=`gmt gmtinfo -i2 -C track$j.data | gawk '{print $1"/"$2}'`"/0/400"
    rm track$j.profile
#    for i in $s1
#    do
#        r=`echo "scale=3;(6378137-$i*1000)/1000"|bc`
#        gmt grdtrack track$j.data -G$i.nc | gawk -v PREC=100 -v dep=$r '{print $2, dep, $4}' >> track_deg$j.profile
#    done
    for i in $s0
    do
        gmt grdtrack track$j.data -G$i.nc | gawk -v PREC=100 -v dep=$i '{print $3, dep, $4}' >> track$j.profile
    done
    
    gmt surface track$j.profile -R${R_interp_cartesian[j]} -I70/36 -V -Gslice$j
    #echo "gmt surface track$j.profile -R${R_interp_cartesian[j]} -I100/50 -V -Gslice$j"

    #extract elevation
    gmt grdtrack track_e$j.data -G$topo >track_e$j.profile
    
    #tac命令将文件逆序输出
    tac track_e$j.data | gawk '{print $1,$2,$3,0}'>>track_e$j.profile


    gawk '{printf"%15.6f %15.6f\n",$3,$4/1000.0}' track_e$j.profile  > elevation$j    


    #Extract depth to the moho along profiles, with respect to sea level
    gmt grdtrack track$j.data -Gdepthtomoho.nc > moho$j.profile
    tac moho$j.profile | gawk '{print $1,$2,$3,0}'>>moho$j.profile
    gawk '{printf"%15.6f %15.6f\n",$3,-1.0*$4}' moho$j.profile  > moho$j
done
for i in $s1
do
    rm $i.nc
done

gmt begin vslice_ini eps E300
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 1p,black
	#gmt set MAP_FRAME_TYPE fancy+
#	gmt makecpt -Cwysiwyg -T-90/90/10 -D -Ic -H > mycpt.cpt
	gmt makecpt -Cwysiwyg -T-60/60/0.1 -Z  -Ic -D -H > mycpt.cpt	
#	gmt grd2cpt slice0 slice1 slice2 -Ic -D -H -Cwysiwyg -E30 >mycpt.cpt
    
    #A
	gmt basemap -Jx0.005/-0.005 -R${img_Region_cartesian[0]} -Bxaf -Bya200f100g100+l"Depth (km)" -BWtsr
	gmt grdimage slice0  -E300 -Q -Cmycpt.cpt
	gmt plot moho0 -W0.5p,black,dashed  -Gwhite -t30
	gmt plot moho0 -Wthick,black,dashed
	#gmt colorbar -DJBC+h+w5c/0.18c+o0.c/1.3c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt
	
	Ra=(`gmt gmtinfo -I- -C elevation0`)
	ele0=${Ra[2]}
	ele1=`echo "${Ra[3]}+2.2"|bc`
	gmt basemap -Jx0.005/0.06 -R${Ra[0]}/${Ra[1]}/$ele0/$ele1 -BWrb -Bya6f2 -Y2c
	gmt plot elevation0 -W0.5p,black  -Ggray

	
	gmt basemap -Jx0.005/0.08 -R${Ra[0]}/${Ra[1]}/$ele0/`echo "$ele1+2"|bc` -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 A
1350 7.6 9p 0 A\234
EOF
	echo 451.82 4.8 -90 0.32 | gmt plot -Sv0.14c+eA -W0.6p -Gblack
	echo 451.82 6.55 8p 0 MFT | gmt text -F+f+A
	echo 800 8.3 8p 0 Tibetan plateau | gmt text -F+f+A
	echo 1250 3.2 8p 0 Tarim | gmt text -F+f+A
	
	
	#B
	Rb=(`gmt gmtinfo -I- -C elevation1`)
	gmt basemap -Jx0.005/-0.005 -R${img_Region_cartesian[1]} -Bxaf -BWtsr -X7.2c -Y-2c
	gmt grdimage slice1  -E300 -Q -Cmycpt.cpt
	gmt plot moho1 -W0.5p,black,dashed  -Gwhite -t30
	gmt plot moho1 -Wthick,black,dashed
	
	gmt basemap -Jx0.005/0.06 -R${Rb[0]}/${Rb[1]}/$ele0/$ele1 -Blrb    -Y2c
	gmt plot elevation1 -W0.5p,black  -Ggray
	gmt basemap -Jx0.005/0.08 -R${Rb[0]}/${Rb[1]}/$ele0/`echo "$ele1+2"|bc` -Bb
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
	Rc=(`gmt gmtinfo -I- -C elevation2`)
	gmt basemap -Jx0.005/-0.005 -R${img_Region_cartesian[2]} -Bxaf+l"Distance (km)" -Bya200f100g100+l"Depth (km)" -BWtSr -X-7.2c -Y-5.25c
	gmt grdimage slice2  -E300 -Q -Cmycpt.cpt
	gmt plot moho2 -W0.5p,black,dashed  -Gwhite -t30
	gmt plot moho2 -Wthick,black,dashed
	
	gmt basemap -Jx0.005/0.06 -R${Rc[0]}/${Rc[1]}/$ele0/$ele1 -BWrb -Bya6f2  -Y2c
	gmt plot elevation2 -W0.5p,black  -Ggray
	gmt basemap -Jx0.005/0.08 -R${Rc[0]}/${Rc[1]}/$ele0/`echo "$ele1+2"|bc` -Bb
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
	gmt basemap -Jx0.005/-0.005 -R${img_Region_cartesian[3]} -Bxaf+l"Distance (km)" -BWtSr -X7.2c -Y-2c
	gmt grdimage slice3  -E300 -Q -Cmycpt.cpt
	gmt plot moho3 -W0.5p,black,dashed  -Gwhite -t30
	gmt plot moho3 -Wthick,black,dashed
	gmt colorbar -Dn-0.38/-0.8+h+w5c/0.2c+o0.c/0c+e -Bxaf -By+L"kg/m@+3@+" -Cmycpt.cpt	
	
	Rd=(`gmt gmtinfo -I- -C elevation3`)
	gmt basemap -Jx0.005/0.06 -R${Rd[0]}/${Rd[1]}/$ele0/$ele1 -BWrb  -Y2c
	gmt plot elevation3 -W0.5p,black  -Ggray
	gmt basemap -Jx0.005/0.08 -R${Rd[0]}/${Rd[1]}/$ele0/`echo "$ele1+2"|bc` -Bb
gmt text -F+f+A <<"EOF"
50 7.6 9p 0 D
1350 7.6 9p 0 D\234
EOF
	echo 596.20 5.0 -90 0.35| gmt plot -Sv0.15c+eA -W0.6p -Gblack
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
rm track$j.profile
done
