model=../gr_inv_pet.nc
extract_slice_tool=Extract_slice
topo=./tibet.nc
if [ ! -f "$topo" ]; then
    echo "Downloading topography data"
    gmt grdcut @earth_relief_15s.grd -R64/112/19/49 -G$topo
else
    echo "${topo} exists"
fi


topography_Region=65/111/20/48
gmt begin topography jpg E500
gmt set FONT_TAG 9p
gmt set FONT_TITLE 9p
gmt set MAP_TITLE_OFFSET 9p
gmt set FONT 9p
gmt set MAP_FRAME_TYPE plain
gmt set MAP_FRAME_PEN 1p,black
gmt set FORMAT_GEO_MAP ddd:mm:ssF
gmt basemap -JQ88/34/10c -R$topography_Region -BWSEN -Bxa10f5 -Bya5f5
gmt grdimage $topo -Cglobe #-I+a45+nt0.5 #-t10 # -I+d

# ======legend图例
gmt legend -F+pthin,black+gWHITE -DjTR+w5.5c+o0.c/0.c --FONT_ANNOT_PRIMARY=6.5p <<EOF
N 2.5c 3.c
S 0.15c f+l+t 0.3 RED thinner,RED, 0.42c Main Frontal Thrust
S 0.15c c 0.13 MAGENTA thinner,black, 0.38c Potassic magmatism 
G 0.03c
S 0.15c - 0.35 - thinner,YELLOW1,- 0.45c Suture zone
S 0.15c - 0.35 - 0p,white, 0.38c (~15-0 Ma)
G 0.03c
S 0.15c - 0.35 - thinner,BLUE 0.42c Faults
S 0.15c c 0.13 ORANGE thinner,black, 0.38c  Ultrapotassic + adakitic  
G 0.03c
S 0.15c - 0.35 - thinner,PURPLE, 0.45c Thrust belt
S 0.15c - 0.35 - 0p,white, 0.38c magmatism (~30-9 Ma)
EOF

cat >annotations <<EOF
88.1 30.7 6.5p,1 -8 Lhasa Block
86 33.5 6.5p,1 -8 Qiangtang Block
85 28.5 6.5p,1 -19 Himalaya Block
95.0 35.0 6.5p,1 -13. Songpan-Ganzi Block
102.55 30.3 6.5p,BLUE -56.5 Xianshuihe Fault
101.75 23.4 6.5p,BLUE -43.5 Red River Fault
77.90 33.5 6.5p,BLUE -55.5 Karakoram Fault
79.5 36.8 6.5p,BLUE -6 Karakax Fault
97.5 36.0 6.5p,BLUE -11 Kunlun Fault
90 39.3 6.5p,BLUE 13 Altyn Tagh Fault
102.5 36.8 6.5p,BLUE -16 Haiyuan Fault
96.75 23.0 6.5p,BLUE -84.5 Sagaing Fault
94.5 37.8 6.5p, -10 Qaidam Basin
83 40. 6.5p, 0 Tarim Basin
78 38.5 6.5p,PURPLE -21 Western Kunlun
98.5 40.2 6.5p,PURPLE -16 Qilian Shan
81 42.8 6.5p, 10 Tien Shan
106.5 31.0 6.5p 0 Sichuan
106. 30.0 6.5p 0 Basin
75.0 26.0 6.5p 0 Indian Plate
91.5 28.6 6.5p,YELLOW1 0 YZS
92.0 31.5 6.5p,YELLOW1 0 BNS
91.9 34.2 6.5p,YELLOW1 -10 JRS
87.6 37.2 6.5p,YELLOW1 -0 AKMS
94 25 6.5p 18 IBR
82 27 6.5p,red -19 MFT
EOF

start=('74.75/27.79' '80.70/25.48' '85.50/23.9' '91.86/21.78' '76.24/28.14' '84.86/26.42' '87.33/25.51')
end=('78.20/40.0' '83.90/37.5' '88.5/36.2' '94.31/34.18' '79.55/39.60' '84.84/39.33' '97.23/39.64')
start_text=('A' 'B' 'C' 'D' 'E' 'F' 'G')
end_text=("A'" "B'" "C'" "D'" "E'" "F'" "G'")
start_text_loc=('74.20/27.95' '80.02/25.71' '84.83/24.05' '92.62/22.01' '77.15/28.16' '85.35/26.56' '88.39/25.51')
end_text_loc=('77.50/40.08' '83.20/37.42' '87.75/36.21' '94.93/34.03' '80.12/39.18' '85.45/39.06' '97.79/39.29')
loop=$(seq 0 6)
for i in $loop; do
    gmt project -C${start[$i]} -E${end[$i]} -G25 -Q -V >track$j.data
    gmt plot -Wthin,GRAY20, track$j.data
	echo ${start_text_loc[$i]}" 6.5p,GRAY15 1 "${start_text[i]} | gawk '{gsub("/"," ");print}' >>annotations
	echo ${end_text_loc[$i]}" 6.5p,GRAY15 1 "${end_text[i]} | gawk '{gsub("/"," ");print}' >>annotations
done

# ============Main frontal thrust
#gmt plot Eurasia_India_PB.dat -W1.0p,red,-
gmt plot MFT.gmt  -Wthinner,red -Sf1.c/0.13c+l+t+o0.2c -Gred

# ============suture
gmt plot sutures.gmt -Wthinner,YELLOW1,-

# ============thrust belt
gmt plot western_kunlun.gmt -Wthinner,PURPLE -Sf2c/0.3c+l+s+o1
gmt plot qilian_shan.gmt -Wthinner,PURPLE #-Sf2c/0.3c+l+s+o1c
#gmt plot TianShan.gmt -Wthinner,PURPLE #-Sf2c/0.3c+l+s+o1c

#=============Fault
gmt plot kunlun_fault.gmt -Wthinner,BLUE, #-Sf2c/0.3c+l+s+o0.8c
gmt plot Altyn_Tagh_Fault.gmt -Wthinner,BLUE, -Sf2c/0.3c+l+s+o1c
gmt plot Karakax.gmt -Wthinner,BLUE, #-Sf2c/0.3c+l+s+o1c
gmt plot Karakoram_Fault.gmt -Wthinner,BLUE, #-Sf2c/0.3c+l+s+o1c
gmt plot Haiyuan_Fault.gmt -Wthinner,BLUE, #-Sf2c/0.3c+l+s+o1c
gmt plot Xianshuihe_Fault.gmt -Wthinner,BLUE, #-Sf2c/0.3c+l+s+o1c
gmt plot Sagaing_Fault.gmt -Wthinner,BLUE, #-Sf2c/0.3c+l+s+o1c
gmt plot Red_river_fault.gmt -Wthinner,BLUE, #-Sf2c/0.3c+l+s+o1c

#=================Rift
#gmt plot Tangra_Yum_Co_rift.gmt -Wthinner,SPRINGGREEN2,

#IBR=Indo-Burma ranges

# plot magmatism
gmt plot -Sc0.1 ./magmatism/Potassic.txt -Wthinnest,black -GMAGENTA
gmt plot -Sc0.1 ./magmatism/Ultrapotasic_adakitc.txt -Wthinnest,black -GORANGE




gmt text -F+f+A annotations

#YZS: Yarlung-Zangbo suture, BNS: Bangong–Nujiang suture, JS: Jinsha suture
#gmt colorbar -Bxa2000f+l"Elevation (m)"
gmt colorbar -DjBL+h+w4.0c/0.2c+o0.25c/0.45c --MAP_FRAME_PEN=faint --FONT=6.5p --MAP_LABEL_OFFSET=-22p --MAP_TICK_LENGTH_PRIMARY=-2.5p --MAP_ANNOT_OFFSET=2p -Bxa2000f+l"Elevation (m)"   
gmt end show
