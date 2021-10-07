#!/bin/bash
# 检查qe,参数设置，输入文件设置
# 输入文件 VASP文件

#/////////////////常变参数\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# 赝势文件位置
upffile="/home/zhangzw/.app/upf"
# 输出文件位置 
out_dir="./"


# ==================check ============================
if  type "pw.x" > /dev/null; then
    echo "">/dev/null
else
    echo "需要安装qe"
    exit 1
fi


if [ ! $# -gt 4 ]; then
    echo "参数不够，需要5个参数，依次是vasp输入文件，k-points,压力"
    exit 1
else
    echo "vasp文件：$1,k-point：$2,$3,$4,压力：$5"
fi

if [ ! -f "$1" ]; then
    echo "输入vasp文件不存在"
    exit 1
fi
#==================获取cpu核心数，超线程=================
let cpu_nums=$(cat /proc/cpuinfo | grep "processor" | wc -l)/2
#==================赝势================================
# for file in $(ls /home/jack/upf)
# do
#     arr[${#arr[@]}]=$file
# done
# 写一个原子质量和赝势对应的文件，多个元素,添加的元素越多越好
cat << EOF > all_ATOMIC_SPECIES
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

if [ -f "${1%.*}-$5.ATOMIC_POSITIONS" ]; then
    rm ${1%.*}-$5.ATOMIC_POSITIONS
fi
let a=9
for i in $(seq 1 ${#atom_arr[@]})
do
    for j in $(seq 1 ${atomNum_arr[$((i-1))]})
    do
        echo ${atom_arr[$((i-1))]} `sed -n ''$a'p' $1` >> ${1%.*}-$5.ATOMIC_POSITIONS;
        ((a=a+1))
    done
done

if [ -f "${1%.*}-$5.ATOMIC_SPECIES" ]; then
    rm ${1%.*}-$5.ATOMIC_SPECIES
fi

# 数组去重，可能需要吧
atom_arr=($(awk -v RS=' ' '!a[$1]++' <<< ${atom_arr[@]}))
for i in $(seq 1 ${#atom_arr[@]})
do
    cat all_ATOMIC_SPECIES |grep "^${atom_arr[$((i-1))]}" >> ${1%.*}-$5.ATOMIC_SPECIES
done

if [ -f "${1%.*}-$5.bfgs.in" ]; then
    rm ${1%.*}-$5.bfgs.in
fi
#==============k-point=========================
declare -i kx=$2
declare -i ky=$3
declare -i kz=$4

#//////////////结构优化输入文件/////////////////
#                可修改成自己的                /
#                                            /
# ////////////////////////////////////////////
cat << EOF >> ${1%.*}-$5.bfgs.in
 &control
    calculation='vc-relax',
    restart_mode='from_scratch',
    outdir='$out_dir' ,
    pseudo_dir ='$upffile',/
    etot_conv_thr=1.0E-5,
    forc_conv_thr=1.0D-4,
 /
&SYSTEM
    ibrav=0, 
    celldm(1)=1.8897268777743552,
    nat=$((nat-9)),
    ntyp=${#atom_arr[@]},
    ecutwfc =70.0,
    occupations='smearing', 
    smearing='methfessel-paxton', 
    degauss=0.02,
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
`sed -n '1,$p' ${1%.*}-$5.ATOMIC_SPECIES`
CELL_PARAMETERS
`sed -n '3,5p' $1`
ATOMIC_POSITIONS (crystal)
`sed -n '1,$p' ${1%.*}-$5.ATOMIC_POSITIONS`
K_POINTS {automatic}
$kx $ky $kz 0 0 0
EOF
#=============删除临时文件====================
if [ -f "${1%.*}-$5.ATOMIC_POSITIONS" ]; then
    rm ${1%.*}-$5.ATOMIC_POSITIONS
fi

if [ -f "all_ATOMIC_SPECIES" ]; then
    rm all_ATOMIC_SPECIES
fi
#===========结构优化，设置断点重启=================

if [ -f "${1%.*}-$5.bfgs.out" ]; then
    grep "JOB DONE." ${1%.*}-$5.bfgs.out
    if [ $? = 0 ]; then
        echo "存在结构优化重启文件，将从此重启"
    else
        echo "存在结构优化重启文件但未完成，重新运行"
        mpirun -np $cpu_nums pw.x <${1%.*}-$5.bfgs.in> ${1%.*}-$5.bfgs.out   
    fi
else
    mpirun -np $cpu_nums pw.x <${1%.*}-$5.bfgs.in> ${1%.*}-$5.bfgs.out
fi


if [ $? -ne 0 ]; then
    exit 1
else
    echo "结构优化成功！"
fi
# while read line
# do
#     if [ "${line:0:15}" = "CELL_PARAMETERS" ]; then
        
#     fi
# done < ${1%.*}-$5.bfgs.in
start=`grep -n "CELL_PARAMETERS" ${1%.*}-$5.bfgs.out |cut -f1 -d:`
let start=`echo $start|awk '{print $NF}'`
echo $start
end=`grep -n "End final coordinates" ${1%.*}-$5.bfgs.out |cut -f1 -d:`
let end=`echo $end|awk '{print $NF}'`
echo $end

if [ -f "${1%.*}-$5.scf.fit.in" ]; then
    rm ${1%.*}-$5.scf.fit.in
fi

cat << EOF >> ${1%.*}-$5.scf.fit.in
 &control
    calculation='scf',
    restart_mode='from_scratch',
    outdir='$out_dir' ,
    pseudo_dir ='$upffile',/
    etot_conv_thr=1.0E-5,
    forc_conv_thr=1.0D-4,
 /
&SYSTEM
    ibrav=0, 
    celldm(1)=1.8897268777743552,
    nat=$((nat-9)),
    ntyp=${#atom_arr[@]},
    ecutwfc =70.0,
    occupations='smearing', 
    smearing='methfessel-paxton', 
    degauss=0.02,
    la2F = .true.,
 /
&ELECTRONS
    conv_thr = 1.0d-12,
    mixing_beta = 0.7,
 /
 &IONS
 /
 ATOMIC_SPECIES
`sed -n '1,$p' ${1%.*}-$5.ATOMIC_SPECIES`
`sed -n ''$start','$((end-1))'p' ${1%.*}-$5.bfgs.out`

K_POINTS {automatic}
$kx $ky $kz 0 0 0
EOF

if [ -f "${1%.*}-$5.scf.fit.out" ]; then
    grep "JOB DONE." ${1%.*}-$5.scf.fit.out
    if [ $? = 0 ]; then
        echo "存在密度自洽重启文件，将从此重启"
    else
        echo "存在密度自洽重启文件但未完成，重新运行"
        mpirun -np $cpu_nums pw.x <${1%.*}-$5.scf.fit.in> ${1%.*}-$5.scf.fit.out    
    fi
else
    mpirun -np $cpu_nums pw.x <${1%.*}-$5.scf.fit.in> ${1%.*}-$5.scf.fit.out 
fi

if [ $? -ne 0 ]; then
    echo "密度自洽计算失败！退出！"
    exit 1
else
    echo "密度自洽计算成功！"
fi

start=`grep -n "CELL_PARAMETERS" ${1%.*}-$5.bfgs.out |cut -f1 -d:`
let start=`echo $start|awk '{print $NF}'`
echo $start
end=`grep -n "End final coordinates" ${1%.*}-$5.bfgs.out |cut -f1 -d:`
let end=`echo $end|awk '{print $NF}'`
echo $end

kx=$kx/2
ky=$ky/2
kz=$kz/2

if [ -f "${1%.*}-$5.scf.in" ]; then
    rm ${1%.*}-$5.scf.in
fi
cat << EOF >> ${1%.*}-$5.scf.in
 &control
    calculation='scf',
    restart_mode='from_scratch',
    outdir='$out_dir' ,
    pseudo_dir ='$upffile',/
    etot_conv_thr=1.0E-5,
    forc_conv_thr=1.0D-4,
 /
&SYSTEM
    ibrav=0, 
    celldm(1)=1.8897268777743552,
    nat=$((nat-9)),
    ntyp=${#atom_arr[@]},
    ecutwfc =70.0,
    occupations='smearing', 
    smearing='methfessel-paxton', 
    degauss=0.02,
    la2F = .true.,
 /
&ELECTRONS
    conv_thr = 1.0d-12,
    mixing_beta = 0.7,
 /
 &IONS
 /
ATOMIC_SPECIES
`sed -n '1,$p' ${1%.*}-$5.ATOMIC_SPECIES`
`sed -n ''$start','$((end-1))'p' ${1%.*}-$5.bfgs.out`

K_POINTS {automatic}
$kx $ky $kz 0 0 0
EOF


if [ -f "${1%.*}-$5.scf.out" ]; then
    grep "JOB DONE." ${1%.*}-$5.scf.out
    if [ $? = 0 ]; then
        echo "存在疏松自洽重启文件，将从此重启"
    else
        echo "存在疏松自洽重启文件但未完成，重新运行"
        mpirun -np $cpu_nums pw.x <${1%.*}-$5.scf.in> ${1%.*}-$5.scf.out        
    fi
else
    mpirun -np $cpu_nums pw.x <${1%.*}-$5.scf.in> ${1%.*}-$5.scf.out 
fi


if [ $? -ne 0 ]; then
    echo "疏松自洽计算失败！退出！"
    exit 1
else
    echo "疏松自洽计算成功！"
fi
kx=$kx/2
ky=$ky/2
kz=$kz/2

if [ -f "${1%.*}-$5.elph.in" ]; then
    rm ${1%.*}-$5.elph.in
fi
cat << EOF >> ${1%.*}-$5.elph.in
&inputph
  tr2_ph=1.0d-12,
  amass(1)= 1.00800001620,
  amass(2)= 32.06000137330,
  amass(3)= 40.07800000000,
  outdir='$out_dir'
  fildvscf='pwscfdv',
  electron_phonon='interpolated',
  alpha_mix=0.2,
  el_ph_sigma=0.0025,
  el_ph_nsigma=20,
  trans=.true.,
  ldisp=.true.,
  nq1=$kx,
  nq2=$ky,
  nq3=$kz,
 /
EOF

if [ -f "${1%.*}-$5.elph.out" ]; then
    grep "JOB DONE." ${1%.*}-$5.elph.out
    if [ $? = 0 ]; then
        echo "存在电声耦合重启文件，将从此重启"
    else
        echo "存在电声耦合重启文件但未完成，重新运行"
        mpirun -np $cpu_nums ph.x <${1%.*}-$5.elph.in> ${1%.*}-$5.elph.out        
    fi
else
    mpirun -np $cpu_nums ph.x <${1%.*}-$5.elph.in> ${1%.*}-$5.elph.out 
fi

if [ $? -ne 0 ]; then
    echo "电声耦合计算失败！退出！"
    exit 1
else
    echo "电声耦合计算成功！"
fi