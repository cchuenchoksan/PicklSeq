minimap2 -a $2 $1 > alignment.sam
samtools view -@ 4 -Sb -o example_alignment.bam alignment.sam
samtools sort -O bam -o sorted_example_alignment.bam example_alignment.bam

samtools bam2fq sorted_example_alignment.bam > minimap_samtools_out.fastq
gzip minimap_samtools_out.fastq
gunzip -c minimap_samtools_out.fastq.gz | chopper -q 1 --minlength 150 --maxlength 1000 | gzip > filtered_reads.fastq.gz
gunzip filtered_reads.fastq.gz

minimap2 -a $2 filtered_reads.fastq > alignment2.sam
samtools view -@ 4 -Sb -o example_alignment2.bam alignment2.sam
samtools sort -O bam -o sorted_example_alignment2.bam example_alignment2.bam
samtools view sorted_example_alignment2.bam > sorted_example_alignment.txt

rm sorted_example_alignment.bam
rm example_alignment.bam
rm alignment.sam

rm minimap_samtools_out.fastq.gz
rm filtered_reads.fastq

rm sorted_example_alignment2.bam
rm example_alignment2.bam
rm alignment2.sam