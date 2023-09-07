#!/bin/bash

temp="$3/temp"

mkdir $temp

echo "$temp/alignment.sam"
minimap2 -a $2 $1 > $temp/alignment.sam
samtools view -@ 4 -Sb -o $temp/temp_alignment.bam $temp/alignment.sam
samtools sort -O bam -o $temp/sorted_temp_alignment.bam $temp/temp_alignment.bam

samtools bam2fq $temp/sorted_temp_alignment.bam > $temp/minimap_samtools_out.fastq
gzip $temp/minimap_samtools_out.fastq
gunzip -c $temp/minimap_samtools_out.fastq.gz | chopper -q 1 --minlength 150 --maxlength 1000 | gzip > $temp/filtered_reads.fastq.gz
gunzip $temp/filtered_reads.fastq.gz

minimap2 -a $2 $temp/filtered_reads.fastq > $temp/alignment2.sam
samtools view -@ 4 -Sb -o $temp/temp_alignment2.bam $temp/alignment2.sam
samtools sort -O bam -o $temp/sorted_temp_alignment2.bam $temp/temp_alignment2.bam
samtools view $temp/sorted_temp_alignment2.bam > sorted_alignment.txt

rm -rf $temp