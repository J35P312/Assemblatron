# Assemblatron - De novo assembly structural variant analysis toolkit

Assemblatron is a toolkit for evaluating De novo assemblies and performing SV calling.

# install

run the install script:
./INSTALL.sh
The install script will download and install fermikit and tiddit. and install the local copy of htslib

Dependencies:
	svdb
	vcftools
	samtools
	python 2.7


# Call SV

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
 
The variants are called using a slighlty modified version of HTSBOX.

# statistics

compute_stats.py
    
    compute various statistics of an assembly.:
        compute_N50.py contigs.bam

    The statistics include N50, L50, assembly size, and the number of contigs.


# De novo assembly and alignment

runFermi.py

    run the fermi assembler and align the contigs to the reference genome:
        python runFermi.py --ref reference.fasta --fastq input.fastq --prefix prefix

contigAln.py

	Align reads to contigs. The resulting bam file may be used for scaffolfing.

# Misc

bam2fq.sh

    convert a bam file to fastq using samtools. the pipeline removes duplicates, supplementary alignments and non-primary alignments.

        ./bam2fq.sh input.bam > output.fastq

