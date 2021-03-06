#!/usr/bin/env bash
# TO BE RUN FROM THE FOLDER CONTAINING THE JSON FILES

sample=$1			# just the numeric ID

zcat S${sample}.json.gz | jq '.samples[0].readGroups[0].metrics[8].x.data[0].values' | tail -n +2 | grep -v ']' | sed 's/,//' > x
zcat S${sample}.json.gz | jq '.samples[0].readGroups[0].metrics[8].y.data[0].values' | tail -n +2 | grep -v ']' | sed 's/,//' > y
size=`cat x | wc -l`
rm -f sample
for id in `seq ${size}`; do
    echo S$sample >> sample
done
paste sample x y
rm sample x y

    
