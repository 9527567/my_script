#!/bin/bash
# 用于moe wash小分子数据库的

# 删除已有的mdb文件
if test -z "$(ls | grep mdb)"; then
	ls
else
	rm *.mdb
fi
which moe > /dev/null 
if [ $? -ne 0 ];then
    echo "请将moe添加至PATH！"
    exit 1
else
    echo "开始执行数据清洗！"
fi
function ceil() {
    floor=$(echo "scale=0;$1/1" | bc -l) 
    add=$(awk -v num1=$floor -v num2=$1 'BEGIN{print(num1<num2)?"1":"0"}')
    echo $(expr $floor + $add)
}

# 获取cpu核心数，超线程
let cpu_nums=$(cat /proc/cpuinfo | grep "processor" | wc -l)/2

if [ $(($cpu_nums % 12)) == 0 ]; then
    split_nums=$(($cpu_nums / 12 ))
else
    split_nums=$(($cpu_nums / 12 + 1))
fi

cat << EOF > split.svl
#svl

function db_Open;
function db_Close;
function db_ImportSD;
function WashMDB;
function db_Minimize;
function db_ExportDB;
function db_Entries;
function db_nEntries;
function split;
function cat;
function tok_cat;
function tok_drop;
function ceil;
function igen;
function length;
function db_Fields;
local function split_mdb[mdb_file,split_nums]
    local entryKey,numsEntryKey,listEntryKey,i,newFile,fields,fields_types;
    entryKey = db_Entries mdb_file;
    numsEntryKey = db_nEntries mdb_file;
    [fields,fields_types] = db_Fields mdb_file;
    listEntryKey = split[entryKey,ceil(numsEntryKey/split_nums)];
    for i in igen split_nums loop 
        newFile = tok_cat[tok_drop[mdb_file,-4],'_',totok(i),'.mdb'];
        db_ExportDB[newFile,mdb_file,fields,listEntryKey(i)];
    endloop
endfunction
local function main[mdb_file,sdf_file,split_nums]
    local fileKey;
    fileKey = db_Open[mdb_file,'create'];
    db_Close fileKey;
    db_ImportSD[mdb_file,sdf_file,'mol'];
    split_mdb[mdb_file,split_nums];
endfunction
EOF
moebatch -exec "run ['split.svl',['tmp.mdb','$1',$split_nums]]"
if test -z "$(ls | grep tmp.mdb)"; then
	ls
else
	rm tmp.mdb
fi
if test -z "$(ls | grep split.svl)"; then
	ls
else
	rm split.svl
fi

cat <<EOF >data_wash.svl
function WashMDB;
function db_Minimize;
local function main[mdb_file]
    local fileKey;
    WashMDB [mdb_file,'mol','',[esel:0,opendbv:0,destfield:'mol',namefrom:'src',namefield:'$File',onecomp:1,salts:1,fragsave:0,fragfield:'salt',neutralize: 0,hydrogens: 'add',protomers : 1,dominantProt: 1,pH: 7 ,original: 0,   enumsize    : 1,enumdup:0,minPctC:0,seqfld:'pseq',depict:0,scale:0,verbose:1]];
    db_Minimize [mdb_file,'mol',[esel:0, dst_field:'', rebuild:'Rebuild3D', filterOpt:[ forceTransAmide:1, forceTransEster:1, forceChair:1, forceTransConjugatedEster:1, forceTransVinyl:1, forceOriginalCC:1, forceOriginalCN:0, forceOriginalNN:0, forceOriginalChirality:0, forceTransUrea:1 ], keep_chirality:0, pot_charge:1, add_h:1, _terse_output:1, gtest:0.1 ]];
endfunction
EOF
for mdb_file in `ls *.mdb`
do 
mdb_array[${#mdb_array[@]}]=${mdb_file}
done

for i in ${mdb_array[@]}
do 
moebatch -exec "run ['data_wash.svl',['$i']]" &
done
wait
if test -z "$(ls | grep data_wash.svl)"; then
	ls
else
	rm data_wash.svl
fi

cat << EOF >combine.svl
# svl
function igen;
function db_ImportDB;
function tok_cat;
function totok;

local function main[mdb_file,new_file,split_nums]
    local i,newName;
    for i in igen split_nums loop
        newName = tok_cat[tok_drop[mdb_file,-4],'_',totok(i),'.mdb'];
        if i == 1 then
            fcopydel[newName,tok_cat[tok_drop[new_file,-4],'.mdb']];
        else
            db_ImportDB[tok_cat[tok_drop[new_file,-4],'.mdb'],newName];
        endif
        fdelete newName;
    endloop
endfunction
EOF
moebatch -exec "run ['combine.svl',['tmp.mdb','$1',$split_nums]]"
if test -z "$(ls | grep combine.svl)"; then
	ls
else
	rm combine.svl
fi