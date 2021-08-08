IFS=''
while read line
do
tmp=`echo $line|awk '{if(NF==5) print $5}'`
#echo $tmp
arr[${#arr[@]}]=$tmp
done < miglitol.molden.chg
while read line
do
tmp=`echo $line|awk '{if(NF==9) print $9}'`
if [ -n "$tmp" ];then
tempval=${arr[i]}
let i=${i}+1
echo "$line"|awk '{if(NF==9) sub($9,""'$tempval'"");print $0}'>>test.mol2
else
echo "$line">>test.mol2
fi
done < miglitol.mol2
