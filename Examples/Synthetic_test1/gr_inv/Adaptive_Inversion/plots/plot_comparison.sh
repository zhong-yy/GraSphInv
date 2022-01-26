model=../ada_result.nc
rslice_adaptive=./test_model_r_slice.nc
lonslice1_adaptive=./test_model_lon_slice1.nc
lonslice2_adaptive=./test_model_lon_slice2.nc

extract_slice_tool=Extract_slice
$extract_slice_tool $model $rslice_adaptive r 6238137
$extract_slice_tool $model $lonslice1_adaptive lon 96.5
$extract_slice_tool $model $lonslice2_adaptive lon 103.5

model2=../../Non-adaptive_Inversion/non_ada_result.nc
rslice_regular=./test_model_r_slice_regular.nc
lonslice1_regular=./test_model_lon_slice1_regular.nc
lonslice2_regular=./test_model_lon_slice2_regular.nc



$extract_slice_tool $model2 $rslice_regular r 6238137
$extract_slice_tool $model2 $lonslice1_regular lon 96.5
$extract_slice_tool $model2 $lonslice2_regular lon 103.5

Region1=90/110/20/40
Region2=90/110/5978.137/6378.137
Region3=20/40/5978.137/6378.137

#gmt begin slices jpg E300
#	gmt set FONT_TAG 15p
#	gmt set FONT_TITLE 12p
	#gmt set FONT 10p
#	gmt set MAP_FRAME_TYPE plain #fancy+
#	gmt subplot begin 2x2 -Fs6.c/4.c -A+JTL+o-0.0c/-0.c -M1.c/1c # -SRl -SCb # -BWSrt

#		gmt subplot set 0 -A'(a)'
#		gmt basemap -JL100/30/25/35/? -R$Region1 -BWSen+t"400km" -Bafg
#		gmt grdimage $rslice_adaptive  -E100 -t10 -Q 
#		gmt colorbar -DJBC+h+w5c/0.2c+o0.c/1.0c -Bxaf -By+L"kg/m@3@"
		
#		gmt subplot set 1 -A'(b)'
#		gmt basemap -JPa?/100z -R$Region2 -Bafg -BWNse+t"30\260 N"
#		gmt grdimage  $latslice  -E100 -t10 -Q 
#		gmt colorbar -DJBC+h+w5c/0.2c+o0.c/1.0c -Bxaf -By+L"kg/m@3@"        	        	        
			        

#		gmt subplot set 3 -A'(c)'
#		gmt basemap -JPa?/30z -R$Region3 -Bafg -BWNse+t"110\260 E"
#		gmt grdimage  $lonslice  -E100 -t10 -Q 
#		gmt colorbar -DJBC+h+w5c/0.2c+o0.c/1.0c -Bxaf -By+L"kg/m@3@"
#	gmt subplot end
	#gmt colorbar -DJBC+h+w7c/0.2c+o0.c/1.0c+e -Bxaf -By+L"kg/m@3@" -Cmycpt.cpt
#gmt end show

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

interpolation=n
gmt begin synthetic_test_gr_model_no_cst eps E300
#	gmt set FONT_TITLE 9p
#	gmt set MAP_TITLE_OFFSET 9p 
	gmt set FONT 9p
	gmt set MAP_FRAME_TYPE plain
	gmt set MAP_FRAME_PEN 0.5p,black
	gmt set MAP_ANNOT_MIN_SPACING 12p
	#gmt set FORMAT_GEO_MAP ddd:mm:ssF
	# gmt set IO_LONLAT_TOGGLE true
        gmt makecpt -Cwysiwyg  -T-90/90 -D -Ic -H > mycpt.cpt #wysiwyg
        mycpt=mycpt.cpt
        
        gmt basemap -JL100/30/25/35/3.8c -R$Region1 -BWSen -Bafg 
		gmt grdimage $rslice_adaptive  -E600 -n$interpolation -Q -C$mycpt #-BWSen -Bafg
		#gmt colorbar -DJMR+w4.c/0.2c+o0.5c/0c -Bxaf -By+L"kg/m@+3@+" -C$mycpt
		
		echo "Depth=140 km" | gmt text -F+cBC+f9p, -D0c/-0.85c -N
		echo "Adaptively refined mesh" | gmt text -F+cTC+f10p, -D0c/0.75c -N
		echo "(a)" | gmt text -F+cTL+f10p, -D-1.5c/0c -N
		gmt plot $true_model_hori -W0.75p,white,
		gmt plot -W1.0p,BROWN3,dashed <<EOF
96.5 20
96.5 40
>
103.5 20
103.5 40
EOF
        echo "95.3 20.8 A" | gmt text -F+f10p,Helvetica-Bold,BROWN3 -N
        echo "95. 39.1 A\234" | gmt text -F+f10p,Helvetica-Bold,BROWN3 -N
        echo "104.6 20.8 B" | gmt text -F+f10p,Helvetica-Bold,BROWN3 -N
        echo "104.9 39.1 B\234" | gmt text -F+f10p,Helvetica-Bold,BROWN3 -N
		
		gmt basemap -JPa4.8c/30z -R$Region3  -BWNse  -Bxafg -Bya200f50g100 -Y-3c -X-0.5c
		gmt grdimage  $lonslice1_adaptive  -E600 -n$interpolation -Q -C$mycpt #-BWNse -Bxafg -Bya200f50g100+l"km"  #nn最邻近插值
		gmt plot $true_model_lon1 -W0.75p,white,
		echo "Longitude=96.5\260 E" | gmt text -F+cBC+f9p, -D0c/-0.36c -N
		echo "(b)" | gmt text -F+cTL+f10p, -D-1.0c/0.35c -N
		echo "A" | gmt text -F+cBL+f10p,Helvetica-Bold,BROWN3 -N -D-0.14c/-0.45c
		echo "A\234" | gmt text -F+cBR+f10p,Helvetica-Bold,BROWN3 -N -D0.18c/-0.45c
		
		gmt basemap -JPa4.8c/30z -R$Region3 -Bxafg -Bya200f50g100 -BWNse -Y-2.6c
		gmt grdimage  $lonslice2_adaptive  -E600 -n$interpolation -Q -C$mycpt #-Bxafg -Bya200f50g100 -BWNse
		#gmt colorbar -DJMR+w4.3c/0.2c+o1c/1.9c+e -Bxa20f10 -By+L"kg/m@+3@+" -C$mycpt
		gmt plot $true_model_lon2 -W0.75p,white,
		echo "Longitude=103.5\260 E" | gmt text -F+cBC+f9p, -D0c/-0.36c -N
		echo "(c)" | gmt text -F+cTL+f10p, -D-1.0c/0.35c -N
		echo "B" | gmt text -F+cBL+f10p,Helvetica-Bold,BROWN3 -N -D-0.15c/-0.45c
		echo "B\234" | gmt text -F+cBR+f10p,Helvetica-Bold,BROWN3 -N -D0.21c/-0.45c		

		
        gmt basemap -JL100/30/25/35/3.8c -R$Region1 -BwSEn -Bafg -X+6.25c -Y+5.6c
		gmt grdimage $rslice_regular  -E600 -n$interpolation -Q -C$mycpt #-BWSen -Bafg
		#gmt colorbar -DJMR+w4.c/0.2c+o0.5c/0c -Bxaf -By+L"kg/m@+3@+" -C$mycpt
		echo "Depth=140 km" | gmt text -F+cBC+f9p, -D0c/-0.85c -N
		echo "Uniform mesh" | gmt text -F+cTC+f10p, -D0c/0.75c -N
		gmt plot $true_model_hori -W0.75p,white,
		gmt plot -W1.0p,BROWN3,dashed <<EOF
96.5 20
96.5 40
>
103.5 20
103.5 40
EOF
        echo "95.3 20.8 A" | gmt text -F+f10p,Helvetica-Bold,BROWN3 -N
        echo "95. 39.1 A\234" | gmt text -F+f10p,Helvetica-Bold,BROWN3 -N
        echo "104.6 20.8 B" | gmt text -F+f10p,Helvetica-Bold,BROWN3 -N
        echo "104.9 39.1 B\234" | gmt text -F+f10p,Helvetica-Bold,BROWN3 -N
        
		gmt basemap -JPa4.8c/30z -R$Region3  -BwNsE  -Bxafg -Bya200f50g100 -Y-3c -X-0.5c
		gmt grdimage  $lonslice1_regular  -E600 -n$interpolation -Q -C$mycpt #-BWNse -Bxafg -Bya200f50g100+l"km"  #nn最邻近插值
		gmt plot $true_model_lon1 -W0.75p,white,
		echo "Longitude=96.5\260 E" | gmt text -F+cBC+f9p, -D0c/-0.36c -N
		echo "A" | gmt text -F+cBL+f10p,Helvetica-Bold,BROWN3 -N -D-0.14c/-0.45c
		echo "A\234" | gmt text -F+cBR+f10p,Helvetica-Bold,BROWN3 -N -D0.18c/-0.45c		
		
		gmt basemap -JPa4.8c/30z -R$Region3 -Bxafg -Bya200f50g100 -BwNsE -Y-2.6c
		gmt grdimage  $lonslice2_regular  -E600 -n$interpolation -Q -C$mycpt #-Bxafg -Bya200f50g100 -BWNse
		gmt colorbar -DJBC+w5.5c/0.25c+o-3.c/0.9c+e -Bxa20f10 -By+L"kg/m@+3@+" -C$mycpt
		gmt plot $true_model_lon2 -W0.75p,white,
		echo "Longitude=103.5\260 E" | gmt text -F+cBC+f9p, -D0c/-0.36c -N
		echo "B" | gmt text -F+cBL+f10p,Helvetica-Bold,BROWN3 -N -D-0.15c/-0.45c
		echo "B\234" | gmt text -F+cBR+f10p,Helvetica-Bold,BROWN3 -N -D0.21c/-0.45c			    
gmt end show
rm $rslice_adaptive
rm $lonslice1_adaptive
rm $lonslice2_adaptive
rm $rslice_regular
rm $lonslice1_regular
rm $lonslice2_regular
rm $true_model_hori
rm $true_model_lon1
rm $true_model_lon2
rm mycpt.cpt


