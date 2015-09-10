#!/usr/bin/env python
# Execute this script where the *.soft GSM metadata files are; requires geo_fields.py

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

cat <(echo -e "GEO accession\tcell type\ttissue\tcell line\tbarcode\tflowcell\tlane\tlibrary id\tplatform\tseqc sample\tsite\ttreatment\tpatient id\tgender\tage\tsubtype\tdisease state") <(for i in *.soft; do paste <(echo $i) <(cat $i | python /home/student/anellor1/SRAmetadata/geo_fields.py); done)