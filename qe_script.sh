#!/bin/bash
# 检查qe,参数设置，输入文件设置
# 输入文件 VASP文件
# ==================check ============================
if  type "pw.x" > /dev/null; then
    echo "">/dev/null
else
    echo "需要安装qe"
    exit 1
fi


if [ ! $# -gt 4 ]; then
    echo "参数不够，需要5个参数"
    exit 1
fi

#==================赝势================================
# for file in $(ls /home/jack/upf)
# do
#     arr[${#arr[@]}]=$file
# done
# 写一个原子质量和赝势对应的文件，多个元素
cat << EOF > ATOMIC_SPECIES
H     1.00800001620       H.pbe-rrkjus_psl.1.0.0.UPF
C    12.01099967960       C.pbe-n-kjpaw_psl.1.0.0.UPF
S    32.06000137330       s_pbe_v1.4.uspp.F.UPF
La   138.90547000000      La.GGA-PBE-paw-v1.0.UPF
Tl   203.973865           Tl_pbe_v1.2.uspp.F.UPF
Ca   40.07800000000       Ca_pbe_v1.uspp.F.UPF
Cu   63.546               Cu_pbe_v1.2.uspp.F.UPF
Sr   87.62                Sr_pbe_v1.uspp.F.UPF
O    15.999               O.pbe-n-kjpaw_psl.0.1.UPF
EOF

#=================bfgs================================

let nat=`awk 'END{print NR}' $1`
for atom in $(sed -n '6p' $1)
do
    atom_arr[${#atom_arr[@]}]=$atom
done
for atomNum in $(sed -n '7p' $1)
do
    atomNum_arr[${#atomNum_arr[@]}]=$atomNum
done

if [ -f "temp.ATOMIC_POSITIONS" ]; then
    rm temp.ATOMIC_POSITIONS
fi
let a=9
for i in $(seq 1 ${#atom_arr[@]})
do
    for j in $(seq 1 ${atomNum_arr[$((i-1))]})
    do
        echo ${atom_arr[$((i-1))]} `sed -n ''$a'p' $1` >> temp.ATOMIC_POSITIONS;
        ((a=a+1))
    done
done
# 数组去重，可能需要吧
atom_arr=($(awk -v RS=' ' '!a[$1]++' <<< ${atom_arr[@]}))
for i in $(seq 1 ${#atom_arr[@]})
do
    cat ATOMIC_SPECIES |grep "^${atom_arr[$((i-1))]}" >> temp.ATOMIC_SPECIES
done

if [ -f "temp.bfgs.in" ]; then
    rm temp.bfgs.in
fi
#//////////////结构优化输入文件/////////////////
#                可修改成自己的               /
#                                           /
# ///////////////////////////////////////////
cat << EOF >> temp.bfgs.in
 &control
    calculation='vc-relax',
    restart_mode='from_scratch',
    outdir='/tmp' ,
    pseudo_dir ='/home/jack/upf',/
    etot_conv_thr=1.0E-5,
    forc_conv_thr=1.0D-4,
 /
&SYSTEM
    ibrav=0, celldm(1)=1.8897268777743552,
    nat=$((nat-9)),
    ntyp=${#atom_arr[@]},
    ecutwfc =70.0,
    occupations='smearing', smearing='methfessel-paxton', degauss=0.02,
    la2F = .true.,
 /
&ELECTRONS
    conv_thr = 1.0d-12,
 /
 &IONS
 /
 &CELL
   cell_dynamics = 'bfgs',
   press = $5,
   cell_factor=1.8,
 /
ATOMIC_SPECIES
`sed -n '1,$p' temp.ATOMIC_SPECIES`
CELL_PARAMETERS
`sed -n '3,5p' $1`
ATOMIC_POSITIONS (crystal)
`sed -n '1,$p' temp.ATOMIC_POSITIONS`
K_POINTS {automatic}
$2 $3 $4 0 0 0
EOF
#=============删除临时文件====================
if [ -f "temp.ATOMIC_POSITIONS" ]; then
    rm temp.ATOMIC_POSITIONS
fi
if [ -f "temp.ATOMIC_SPECIES" ]; then
    rm temp.ATOMIC_SPECIES
fi
if [ -f "ATOMIC_SPECIES" ]; then
    rm ATOMIC_SPECIES
fi