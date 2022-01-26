reference_model="SL2013sv_25k-0.5d.mod"
#sed '1d' $reference_model | gawk -v west=60 -v east=116 -v south=16 -v north=52 -v dep0=0 -v dep1=600 '$1<=dep1 && $1>=dep0 && $2>=west && $2<=east && $3>=south && $3<=north {printf "%20f %20f %20f %20f\n",$3,$2,$1,$6}' $reference_model > Tibetan_dVs
sed '1d' $reference_model | gawk -v west=60 -v east=116 -v south=16 -v north=52 -v dep0=0 -v dep1=600 '$1<=dep1 && $1>=dep0 && $2>=west && $2<=east && $3>=south && $3<=north {printf "%20f %20f %20f %20f\n",$3,$2,$1,$7}' $reference_model > Tibetan_VsAbs
