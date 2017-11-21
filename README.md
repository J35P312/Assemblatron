# Assemblatron - De novo assembly structural variant detection/validation tool

0: install
run the install script:
./INSTALL.sh
The install script will download and install fermikit and tiddit. and install the local copy of htslib

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

#other scripts

compute_coverage.py
    
    compute the coverage across the genome, as well as the perentage of the genome having 0 coverage.

    1:generate a coverage tab file using tiddit:
	./TIDDIT/bin/TIDDIT --cov -b file.bam -z 50 -o prefix

    2: compute the coverage
	python compute_coverage.py prefix.tab

compute_N50.py
    
    compute the N50 of an assembly:
        compute_N50.py input.fa

contigSupport.py
	accepts contigs in fasta or bam format(suffix, bam,fasta, or fa), as well as reads in fastq or bam format(suffix: bam, fastq, or fq). The reads are mapped to the contigs, and three files are produced (based on --prefix)
	
	prefix.bam - a bam file with the results of mapping the reads to the contigs
	prefix.tab - the coverage across bins on the contigs
	prefix.stats - per contig statistics and dataset statistics


	
 
runFermi.py
    run the fermi assembler and align the contigs to the reference genome:
        python runFermi.py --ref reference.fasta --fastq input.fastq --prefix prefix
    fermikit is assumed to be installed in the Assemblatron folder, if this is not the case, you need to change the fermi variable in the runFermi script
    alternatively, the pipeline can

    runFermi.py will not apply any error correction of the reads, nor any trimming. By skipping these steps, the run time is greatly reduced. Additionally, the settings are less stringent compared to standard fermikit.
    runfermi.py will not perform any variant calling: the end result is a bam file.
 
bam2fq.sh
    convert a bam file to fastq using samtools. the pipeline removes duplicates, supplementary alignments and non-primary alignments.

        ./bam2fq.sh input.bam > output.fastq

    
#De novo assembly

Assemblatron has been tested using 10X supernova output, the "raw" fasta file has been found to give the best results (high sensitivity).
The aligned FermiKit  output works well with Assemblatron aswell. Such file can be created using FermiKit or the runFermi.py wrapper

