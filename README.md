# Assemblatron - De novo assembly tool kit

Assemblatron is a toolkit for De novo assembly, scaffolding, and variant calling of de  novo assemblies.  type 
	python assemblatron.py --help

For a help message of the avialable tools.

# The workflow 

The full workflow involves:

    1: assembly using Fermi2 (similar to workflow as fermiKit)

    2: scaffolding using BESST (using the same input data as for 1, or a mate pair library)

    3: alignment of the contigs

    4: Quality control, using assemblatron stats or quast

    5: variant calling

However, you may start or end the analysis as any step.
The assemblatron workflow is run through the following commands:

# Install

run the install script:
./INSTALL.sh
The install script will download and install BESST, fermikit,kmc, bfc-kmc, and quast.
Dependencies:

	conda
	samtools
	python 2.7
	bwa
	numpy

Assemblatron may also be run using the singularity container:



# Assemble
Assemblatron performs de novo assembly using a pipeline similar to fermikit. The pipeline is run through the following commands:

		python assemblatron.py --assemble --fastq <input.fastq> --prefix <prefix>

You may also align the contigs to the reference genome:

		python assemblatron.py --assemble --ref <reference.fasta> --fastq <input.fastq> --prefix <prefix> --align

Per default, the Assembly module will generate an output file named <prefix>.fastq. If asssemblatron was run with the --align option, a bam file named <prefix>.bam will be produced as well.

Type help for more  info:

		 python assemblatron.py --assemble --help

The main differences between fermikit and Assemblatron is that Assemblatron uses the kmc branch of bfc. Additionally, Assemblatron uses piping more extensively, hence it  will not create any intermidiate fastq files.
For more info on Fermikit, visit the Fermikit website:

		https://github.com/lh3/fermikit

# Scaffolfing

Perform Scaffolding using BESST. 

		python assemblatron.py scaffold --contigs <contigs> --fastq <fastq> --output <output> --tmpdir
	
Where <Contigs> is a fasta file containing the contigs produced through de novo assembly, and <fastq> is a fastq file containing paired-end reads. The scaffolding process including alignment takes about 24-48 hours on a 16 core machine when using a 30X whole genome for scaffolding.
When running BESST on the assemblatron output contigs, the N50 is usually improved by a factor of 10 or so.

# Align
Assemblatron performs alignment using bwa mem. The output is printed to a file named <prefix.bam>. The input is asssumed to be contigs produced through any assembly process.

	python assemblatron.py --align --ref <reference.fasta> --contigs <contigs.fasta> --prefix <prefix>

A bam file named <prefix>.bam will be produced

# Stats
compute various statistics of an assembly. The output is printed to stdout.

python assemblatron.py --stats --bam <contigs_bam>
	
The statistics include N50, L50, assembly size, and the number of contigs.

# Quast

The quality control may be  performed using quast. Quast performs a more in-depth but slower analysis.

python assemblatron.py --quast --fasta <contigs_fasta> --output <output_folder>

python assemblatron.py --quast --fasta <contigs_fasta> --ref <reference.fasta> --output <output_folder>

The statistics module supports any number of fasta files. Type --help for more information. 
NOTE: use absolute path for the output directory.

# SV
Call SV. Assemblatron will calssify variants as DEL, INV, BND (translocation or complex), INS, and TDUP.

    python assemblatron.py --sv --bam <contigs.bam> > out.vcf

The output is  printed to stdout

other options:

    -h, --help            show this help message and exit
    --sv                  call SV from the aligned contigs
    --bam BAM             input bam (contigs
    --q Q                 minimum allowed mapping quality(default = 10)
    --len_ctg LEN_CTG     minimum uniquely mapped contig length (default = 40)
    --max_coverage MAX_COVERAGE
    --sample              sample id (default = bam-filename)
                        calls from regions exceeding the maximum coverage are
                        filtered
    --min_size MIN_SIZE   minimum variant size (default=100)

# SNV
Call indels and SNV using htsbox pileup (same as fermikit).

the  output is printed to stdout

options:

  -h, --help  show this help message and exit
  --bam BAM   input bam (contigs)
  --ref REF   reference fasta

# fastq

convert a bam file to fastq.

	python assemblatron.py --fastq <input.bam>

The output is  printed to stdout.

# fasta
Convert a contig bam file to fasta.

	python assemblatron.py --fastq <contigs.bam>

The output is  printed to stdout.

# Testing
De novo assembly of NA12878 may be  downloaded from the 10X website:

Start by aligning the contigs using the align module, thereafter you may call SV and SNV using the sv and snv modules. Truthsets are available at the GIAB  website:

As well as through the SV classify supplementary methods:

# Cite
Cite the components that you used, as well as the Assemblatron git hub page.
 
	If you performed De novo assembly, cite the FermiKit paper:

		https://github.com/lh3/fermikit

	If you performed scaffolding, please cite the BESST paper:

		https://github.com/ksahlin/BESST

	For more info on QUAST:

		https://github.com/ablab/quast
