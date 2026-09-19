#!/bin/bash

if [ ! -f treat.tsv ]; then
    echo "treat.tsv 文件不存在，退出脚本。"
    exit 1
fi

awk 'BEGIN{FS=OFS="\t"}NR==1{print "SampleID","SampleName","Treat"}NR>1{print $5,$5,$4}' treat.tsv > map.tsv
awk 'BEGIN{FS=OFS="\t"}{print $1,$5}' treat.tsv > rename.tsv 

sed "1d" treat.tsv | \
    awk 'BEGIN{FS=OFS="\t"}{ print $5,
        "$PWD/01_rawdata/"$5".R1.fq.gz", 
        "$PWD/01_rawdata/"$5".R2.fq.gz"}' | 
    sed "1iSampleID\tforward-absolute-filepath\treverse-absolute-filepath" \
    > manifest.tsv
