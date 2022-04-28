rslice=./test_model_r_slice.nc
lonslice1=./test_model_lon_slice1.nc
lonslice2=./test_model_lon_slice2.nc
xyinc1=81+n/41+n
xyinc2=81+n/81+n
Region1=90/110/20/40
Region3=20/40/5978.137/6378.137

python getSlice.py
gmt xyz2grd "./interp_data1.xyz" -I$xyinc1 -R$Region3 -G$lonslice1
gmt xyz2grd "./interp_data2.xyz" -I$xyinc1 -R$Region3 -G$lonslice2
gmt xyz2grd "./interp_data3.xyz" -I$xyinc2 -R$Region1 -G$rslice

block_dep=(6178.137 6298.137)
block_lat=(27 33)
block_lon=(95 98)

block_dep2=(6178.137 6298.137)
block_lat2=(32 35)
block_lon2=(102 105)

block_dep3=(6178.137 6298.137)
block_lat3=(25 28)
block_lon3=(102 105)


true_model_hori='true_mode_hori'
true_model_lon1='true_model_lon1'
true_model_lon2='true_model_lon2'

if [ -f "./${true_model_hori}" ];then
echo "delete ${true_model_hori}"
rm $true_model_hori
fi

if [ -f "./${true_model_lon1}" ];then
echo "delete ${true_model_lon1}"
rm $true_model_lon1
fi

if [ -f "./${true_model_lon2}" ];then
echo "delete ${true_model_lon2}"
rm $true_model_lon2
fi


python make_frame.py ${block_lon[0]} ${block_lon[1]} ${block_lat[0]} ${block_lat[1]} $true_model_hori
python make_frame.py ${block_lon2[0]} ${block_lon2[1]} ${block_lat2[0]} ${block_lat2[1]} $true_model_hori
python make_frame.py ${block_lon3[0]} ${block_lon3[1]} ${block_lat3[0]} ${block_lat3[1]} $true_model_hori

python make_frame.py ${block_lat[0]} ${block_lat[1]} ${block_dep[0]} ${block_dep[1]} $true_model_lon1


python make_frame.py ${block_lat2[0]} ${block_lat2[1]} ${block_dep2[0]} ${block_dep2[1]} $true_model_lon2
python make_frame.py ${block_lat3[0]} ${block_lat3[1]} ${block_dep3[0]} ${block_dep3[1]} $true_model_lon2

gmt begin synthetic_test_gr_model_no_cst jpg E300
#	gmt set FONT_TITLE 9p
#	gmt set MAP_TITLE_OFFSET 9p 
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 0.25p,black
	#gmt set FORMAT_GEO_MAP ddd:mm:ssF
	# gmt set IO_LONLAT_TOGGLE true
        gmt makecpt -Cwysiwyg  -T-85/85  -D -Ic -H > mycpt.cpt #Cwysiwyg
        mycpt=mycpt.cpt
		gmt basemap -JPa4.9c/30z -R$Region3 -Bxafg -Bya200f50g100 -BWNse
		gmt grdimage  $lonslice2  -E600 -nc -Q -C$mycpt -Bxafg -Bya200f50g100 -BWNse
		gmt colorbar -DJMR+w4.3c/0.2c+o1c/1.9c+e -Bxa20f10 -By+L"kg/m@+3@+" -C$mycpt
		gmt plot $true_model_lon2 -W1.25p,white,
		echo "(c) Longitude=103.5\260 E" | gmt text -F+cTL+f9p, -D0.55c/1.2c -N        	        
			        
		gmt basemap -JPa4.9c/30z -R$Region3  -BWNse -Y+3.5c -Bxafg -Bya200f50g100 #+l"km"
		gmt grdimage  $lonslice1  -E600 -nc -Q -C$mycpt -BWNse -Bxafg -Bya200f50g100 #+l"km"
		gmt plot $true_model_lon1 -W1.25p,white,
		echo "(b) Longitude=96.5\260 E" | gmt text -F+cTL+f9p, -D0.55c/1.2c -N  		

	
		gmt basemap -JL100/30/25/35/3.8c -R$Region1 -BWSen -Bafg -Y-3c -X-5.4c
		gmt grdimage $rslice  -E600 -nc -Q -C$mycpt -BWSen -Bafg
		#gmt colorbar -DJMR+w4.c/0.2c+o0.5c/0c -Bxaf -By+L"kg/m@+3@+" -C$mycpt
		echo "(a) Depth=140 km" | gmt text -F+cTL+f9p, -D0.5c/0.9c -N
		gmt plot $true_model_hori -W1.25p,white,
gmt end show

rm $rslice
rm $lonslice1
rm $lonslice2

rm $true_model_hori
rm $true_model_lon1
rm $true_model_lon2

rm interp_data1.xyz
rm interp_data2.xyz
rm interp_data3.xyz
rm gmt.history
rm mycpt.cpt

