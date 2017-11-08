# Assemblatron - De novo assembly structural variant detection/validation tool

0: install
enter the htsbox directory, and run make
cd htsbox
cd make

Dependencies:
	svdb
	vcftools
	samtools
	python 2.7

1: align the contigs using bwa mem:
        bwa mem -x intractg  reference.fasta contigs.fa | samtools view -bhS - > contigs.bam

2: perform variant calling
        
    python assemblatron.py --bam contigs.bam --ref reference.fasta > calls.vcf

The output file is a vcf, variants are returned as inversion or breakend variant. Assemblatron detects variants of any size and type.
optional parameters:

  --working_dir WORKING_DIR temporary analysis files will be stored here(default=work)
  --bam BAM             input bam
  --q Q                 minimum allowed mapping quality(default = 10)
  --len LEN             minimum length of alignments(default = 200)
  --len_ctg LEN_CTG     minimum contig length(default = 1000)
  --max_dist MAX_DIST   maximum distance of two breakpoints within a contig(default = 100)
 
note: assemblatron is a prototype software.

#Algorithm
A wrapper around a slighly modified version of htsbox.

#Future additions
A variant imputation function will be added
