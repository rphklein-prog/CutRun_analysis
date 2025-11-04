#!/bin/bash

# Cut&Run analysis script

# assumes files are in .fq.gz format, and F and R sequences are located together
# in a subfolder within the data directory. Takes four arguments: the first is the
# directory where the raw data is stored, the second is the directory for bam files,
# the third the location of the bowtie index, the fourth is the chrom sizes file 
# (can be downloaded from UCSC genome browser). 
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <data_directory> <bam_output_dir> <bowtie_index> <chrom_sizes>"
    exit 1
fi

data_directory=$1
bam_file_directory=$2
bowtie_index=$3
chrom_sizes=$4


# adapter sequences can be modified if necessary

adapt_seqs="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT,GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG"

# if bam_file_directory doesn't exist yet, make it

mkdir -p "$bam_file_directory"

# loop through paired files, trim adapters from reads using bbduk, and align using bowtie2 with 
# recommended Cut&Run parameters, output into bam_file_directory

for in1 in "$data_directory"/*/*_1.fq.gz; do
    # infer paths
    in2="${in1/_1.fq.gz/_2.fq.gz}"
    sample_name=$(basename "$(dirname "$in1")")
    echo "Processing sample: $sample_name"
	bbduk.sh in1="$in1" in2="$in2" \
	out1="${sample_name}_clean_1.fq" out2="${sample_name}_clean_2.fq" \
	literal=$adapt_seqs \
	ktrim=r k=23 mink=11 hdist=1 tpe tbo
	
	bowtie2 --local --very-sensitive --no-mixed --no-discordant \
	--phred33 -I 10 -X 700 \
	-x "${bowtie_index}" \
	-1 "${sample_name}_clean_1.fq" -2 "${sample_name}_clean_2.fq" \
	-S $i\.sam
	samtools view -bS "${sample_name}.sam" > "$bam_file_directory/${sample_name}.bam"
done

# print out stats about alignment
for file in $bam_file_directory/*.bam
do
	echo $file
	samtools view -c $file
	samtools view -c -F 260 $file
done

# move to bam_file_directory, for each file, convert to bed, format reads

cd "$bam_file_directory"
for bam in *.bam; do
    base="${bam%.bam}"
    bedtools bamtobed -i "$bam" -bedpe > "${base}.bed"
    awk '$1==$4 && $6-$2 < 1000 {print $0}' "${base}.bed" > "${base}.clean.bed"
    sed -i '/-1/d' "${base}.clean.bed"
    bedSort "${base}.clean.bed" "${base}.clean.sorted.bed"
    bedtools genomecov -bg -i "${base}.clean.sorted.bed" -g "$chrom_sizes" > "${base}.bdg"
done


# to make bigwig files (if desired, otherwise can be commented out)	
for i in `ls -1 *.bdg | sed 's/.bdg//'`
do
	bedSort "$i.bdg" "${i}_sorted.bdg"
	bedClip "${i}_sorted.bdg" "$chrom_sizes" "${i}_sorted_clipped.bdg"
	bedGraphToBigWig "${i}_sorted_clipped.bdg" "$chrom_sizes" "${i}.bw"
done

#for i in `ls -1 *.bam | sed 's/.bam//'`
#do
#echo "$i.bam"
#samtools view $i.bam > $i_sorted.bam
#samtools index $i_sorted.bam
#done