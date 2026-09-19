#!/bin/bash

sed "1d" rename.tsv | while IFS=$'\t' read -r old_name new_name; do
    r1="01_rawdata/${old_name}_R1.fq"
    r2="01_rawdata/${old_name}_R2.fq"
    if [ -e "$r1" ] && [ -e "$r2" ]; then
        mv "$r1" "01_rawdata/${new_name}_R1.fq"
        mv "$r2" "01_rawdata/${new_name}_R2.fq"
        echo "renamed $old_name to $new_name"
    else
        echo "file(s) for $old_name do not exist."
    fi
done