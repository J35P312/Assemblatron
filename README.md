# Assemblatron - De novo assembly tool kit

Assemblatron is a toolkit for De novo assembly, scaffolding, and variant calling of de  novo assemblies.  type 
	python assemblatron.py --help

For a help message of the avialable tools. The full workflow involves:

1: assembly using Fermi2

2: scaffolding using BESST (using the same input data as for 1, or a mate pair library)

3: alignment of the contigs

4: compute assembly statistics

5: variant calling

# Install

run the install script:
./INSTALL.sh
The install script will download and install BESST, fermikit, and tiddit. and install the local copy of htslib

Dependencies:
	svdb
	vcftools
	samtools
	python 2.7
	bwa

# Assemble
Assemblatron performs de novo assembly using a pipeline similar fermikit fermikit. The pipeline is run through the following commands:

		python assemblatron.py --assemble --fastq <input.fastq> --prefix <prefix>

You may also align the contigs to the reference genome:

		python assemblatron.py --assemble --ref <reference.fasta> --fastq <input.fastq> --prefix <prefix> --align

Type help for more  info:

		 python assemblatron.py --assemble --help

The main difference between Fermikit and Assemblatron assemble is the default parameters. Additionally, Assemblatron uses piping more  extensively to reduce the writing/reading to disk.
For more info on Fermikit, visit the Fermikit website:

		https://github.com/lh3/fermikit

# Scaffolfing

Perform Scaffolding using BESST. 

		python assemblatron.py scaffold --contigs <contigs> --fastq <fastq> --output <output> --tmpdir
	
Where <Contigs> is a fasta file containing the contigs produced through de novo assembly, and <fastq> is a fastq file containing paired-end reads. The scaffolding process including alignment takes about 24-48 hours on a 16 core machine when using a 30X whole genome for scaffolding.
When running BESST on the assemblatron output contigs, the N50 is usually improved by a factor of 10 or so.

# Align

# Stats
compute various statistics of an assembly.

python assemblatron.py --stats <contigs_bam>
	
The statistics include N50, L50, assembly size, and the number of contigs.

# SV

# SNV

# Conversion

# Cite
Cite the components that you used, as well as the Assemblatron git hub page.
 
	If you performed De novo assembly, cite the FermiKit paper:

		https://github.com/lh3/fermikit

	If you performed scaffolding, please cite the BESST paper:

		https://github.com/ksahlin/BESST
