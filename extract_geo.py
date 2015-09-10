#!/usr/bin/env python
# Execute this script where the *.soft GSM metadata files are

for i in *.soft; do grep "\!Sample_characteristics_ch1" $i | awk -F ':' '{print $1}' | awk -F '=' '{print $2}' | sed -e 's/^[ \t]*//' | sort | uniq; done | sort | uniq -c | awk '$1 > 400 {print $0}' | >/home/student/anellor1/geq400Fields

# From the results, looks like the most important fields are 
# (sample count) (field name)
#   3795 cell type
#   3180 tissue
#   2911 cell line
#   1845 barcode
#   1767 flowcell
#   1767 lane
#   1720 library id
#   1720 platform
#   1720 seqc sample
#   1720 site
#   1098 treatment
#    749 patient id
#    699 gender
#    691 age
#    499 subtype
#    467 disease state

for i in *.soft; do grep "\!Sample_characteristics_ch1" $i | awk -F '=' -v l=$i '{print $2}' | awk -F ':' '{for(i=1;i<+16;i++){a[i]="NA"} gsub(/^[ \t]+/,"",$1);gsub(/[ \t]+$/,"",$1); if ($1 == "cell type") {a[1]=$2} if ($1 == "tissue") {a[2]=$2} if ($1 == "cell line") {a[3]=$2} if ($1 == "barcode") {a[4]=$2} if ($1 == "flowcell") {a[5]=$2} if ($1 == "lane") {a[6]=$2} if ($1 == "library id") {a[7]=$2} if ($1 == "platform") {a[8]=$2} if ($1 == "seqc sample") {a[9]=$2} if ($1 == "site) {a[10]=$2} if ($1 == "treatment") {a[11]=$2} if ($1 == "patient id") {a[12]=$2} if ($1 == "gender") {a[13]=$2} if ($1 == "age") {a[14]=$2} if ($1 == "subtype") {a[15]=$2} if ($1 == "disease state") {a[16]=$2} printf l; for(i=1;i<=16;i++) {printf "\t" a[i]} printf "\n"; }}'; done