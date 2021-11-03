#!/bin/bash
# 使用：resp_orca.sh example.mol2 result.mol2
# 使用sob老师推荐的泛函基组
keyword="! B3LYP D3 def2-TZVP def2/J RIJCOSX noautostart opt miniprint nopop CPCM(water)"
maxcore="%maxcore 1000"
pal="%pal nprocs 6 end"
ORCA=/home/jack/.app/orca

export inname=$1
filename=${inname%.*}

Multiwfn $1 > /dev/null << EOF
100
2
12
tmp.inp
1
0
q
EOF
sed '1,3d' -i tmp.inp
sed '1i '"$keyword"'' -i tmp.inp
sed '2i '"$maxcore"'' -i tmp.inp
sed '3i '"$pal"'' -i tmp.inp

cat tmp.inp
echo Running optimization task via ORCA..
$ORCA/orca tmp.inp 

echo Running orca_2mkl...
$ORCA/orca_2mkl tmp -molden

echo Running Multiwfn...
Multiwfn tmp.molden.input > /dev/null << EOF
7
18
1
y
0
0
100
2
1
opt.pqr
0
q
EOF

rm tmp.cpcm tmp.densities tmp.inp tmp.gbw tmp.molden.input tmp_property.txt

obabel opt.pqr -ipqr -omol2 -O $2
rm opt.pqr
# IFS=''
# while read line
# do
# tmp=`echo $line|awk '{if(NF==5) print $5}'`
# #echo $tmp
# arr[${#arr[@]}]=$tmp
# done < tmp.molden.chg
# while read line
# do
# tmp=`echo $line|awk '{if(NF==9) print $9}'`
# if [ -n "$tmp" ];then
# tempval=${arr[i]}
# let i=${i}+1
# echo "$line"|awk '{if(NF==9) sub($9,""'$tempval'"");print $0}'>>$2
# else
# echo "$line">>$2
# fi
# done < tmp.mol2

