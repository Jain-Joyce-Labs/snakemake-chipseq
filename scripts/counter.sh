#!/bin/sh

date > $4

for var in $1 $2 $3
do
    echo $var >> $4
    samtools view -c -F 260 $var >> $4
done
