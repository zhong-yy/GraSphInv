function check_and_delete(){
    if [ -f "$1" ];then
        echo "delete $1"
        rm $1
    fi
}


PATH_TO_LITHO10_BIN=./bin

function LAB(){
    while read line; do
        latitude=`echo $line | gawk '{print $2}'`
        longitude=`echo $line | gawk '{print $1}'`
        distance=`echo $line | gawk '{print $3}'`
        array=(${string//,/ })  
    #    echo $coordinates
        LAB=`${PATH_TO_LITHO10_BIN}/access_litho -p $latitude $longitude  | gawk '{print $1/1000.0}' | sed -n '1p'`
        echo -e "${longitude}\t${latitude}\t${distance}\t${LAB}">>$2
        #echo -e "\n"
    done < $1
}

loop=$(seq 0 3)
for j in $loop
do
    check_and_delete LAB$j.profile
    LAB track$j.data LAB$j.profile
    gawk '{printf"%15.6f %15.6f\n",$3,$4}' LAB$j.profile  > LAB$j
done


