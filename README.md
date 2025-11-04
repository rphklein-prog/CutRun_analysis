This script takes raw Cut&Run sequence data through all of the steps to prepare it for SEACR enrichment analysis. This includes adapter trimming (with bbmap tool bbduk), 
alignment (with bowtie2, and recommended settings for Cut&Run), data formatting (sorting, converting file formats), and generation of bedgraph files for SEACR analysis, 
and bigwig files (if desired). It requires the following programs: bbmap, bowtie2, samtools, bedtools, bedSort, bedClip, and bedGraphToBigWig.
