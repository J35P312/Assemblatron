import sys
import os
import argparse
import numpy
import math
import itertools

wd=os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, '{}/scripts'.format(wd))
import call
import stats

def assemble(args,wd):
	fermi="{}/fermikit/fermi.kit/".format(wd)
	bfc_path="{}/bfc".format(wd)

	#run kmc
	kmc_prefix="{}.kmc".format(args.prefix)
	kmc="kmc -k{} -m30 {} {} {}".format(args.k,args.fastq,kmc_prefix,args.tmp)
	os.system(kmc)

	#apply bloom filter and build the index
	ropebwt="{}/ropebwt2 -m {} -dNCr - > {}.fmd 2> {}.fmd.log".format(fermi,args.batch,args.prefix,args.prefix)
	bfc="{}/bfc-kmc -k {} -t {} {} {} 2> {}.flt.fq.gz.log".format(bfc_path,args.k,args.cores,kmc_prefix,args.fastq,args.prefix)
	os.system("{} | {}".format(bfc,ropebwt))

	#assemble
	os.system( "{}/fermi2 assemble -l {} -t {} {}.fmd 2> {}.pre.gz.log | gzip -1 > {}.pre.gz".format(fermi,args.l,args.cores,args.prefix,args.prefix,args.prefix) )
	os.system("{}/fermi2 simplify -CS -T {} {}.pre.gz 2>  {}.mag.gz.log > {}.fastq".format(fermi,args.l,args.prefix,args.prefix, args.prefix))

	if args.align:
		os.system("bwa mem -x intractg -t {} {} {}.fastq | samtools view -Sbh - | samtools sort -m 2G - > {}.bam".format(args.cores,args.ref,args.prefix,args.prefix))
		os.system("samtools index {}.bam".format(args.prefix) )

version = "0.2.0"
parser = argparse.ArgumentParser("""Assemblatron: a de novo assembly  pipeline""".format(version),add_help=False)
parser.add_argument('--assemble'       , help="Perform de novo assembly using the Fermi2 assembler", required=False, action="store_true")
parser.add_argument('--scaffold'      , help="perform scaffolding using BESST", required=False, action="store_true")
parser.add_argument('--sv'             , help="call SV from the aligned contigs", required=False, action="store_true")
parser.add_argument('--snv'             , help="call snv from the aligned contigs", required=False, action="store_true")
parser.add_argument('--stats'          , help="compute assembly stats from aligned contigs", required=False, action="store_true")
parser.add_argument('--align'          , help="align contigs to reference using bwa mem", required=False, action="store_true")
parser.add_argument('--fasta'      , help="convert aligned contigs bam file to fasta", required=False, action="store_true")
parser.add_argument('--fastq'      , help="convert bam to fastq", required=False, action="store_true")
parser.add_argument('--quast'        , help="compute assembly stats using quast", required=False, action="store_true")

args, unknown = parser.parse_known_args()

if args.assemble:

	parser = argparse.ArgumentParser("""Assemblatron assemble - a wrapper for the fermi assembler""")
	parser.add_argument('--assemble'       , help="Perform de novo assembly using the Fermi2 assembler", required=False, action="store_true")
	parser.add_argument('--fastq',required = True, type=str, help="input fastq")
	parser.add_argument('--prefix',required = True,type=str, help="prefix of the output files")
	parser.add_argument('--cores',type=int, default =16, help="number of cores (default = 16)")
	parser.add_argument('--batch',type=str, default ="20g", help="batch size for multi-string indexing; 0 for single-string (default=20g)")
	parser.add_argument('-l',type=int, default =81, help="min match (default = 81)")
	parser.add_argument('-k',type=int, default =81, help="minimum kmer length for kmc/bfc error correction (default = 41)")
	parser.add_argument('--align', help="align contigs to reference using bwa mem", required=False, action="store_true")
	parser.add_argument('--ref',type=str, help="reference fasta, required for alignment of the contigs")
        parser.add_argument('--tmp',type=str,default="$TMPDIR", help="tmp directory, kmc will write tmp files here (default=$TMPDIR)")
	args= parser.parse_args()

	if not os.path.isdir(args.tmp):
		print "error: no such directory {}".format(args.tmp)
		print "set the --tmp variable to an existing folder"
		quit()

	if not args.align or (args.align and args.ref):
		assemble(args,wd)
	elif args.align and not args.ref:
		print ("you need a reference to align the contigs: Please supply the reference path through the --ref parameter")

elif args.sv:

	parser = argparse.ArgumentParser("""Assemblatron sv - a variant caller using aligned contigs""")
	parser.add_argument('--sv'        , help="call SV from the aligned contigs", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	parser.add_argument('--q',type=int, default =10 ,help="minimum allowed mapping quality(default = 10)", required=False)
	parser.add_argument('--len_ctg'       ,type=int, default = 40, help="minimum uniqyelly mapped contig length(default = 40)", required=False)
	parser.add_argument('--max_coverage'       ,type=int, default = 5, help="calls from regions exceeding the maximum coverage are filtered", required=False)
	parser.add_argument('--min_size'       ,type=int, default = 100, help="minimum variant size)", required=False)  
        parser.add_argument('--sample'       ,type=str, help="sample id, as shown in the format column", required=False)
	args= parser.parse_args()

	call.main(args)

elif args.stats:
	parser = argparse.ArgumentParser("""Assemblatron stats - compute assembly stats""")
	parser.add_argument('--stats'        , help="compute assembly stats from aligned contigs", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	args= parser.parse_args()

	stats.assembly_stats(args)
elif args.quast:

        parser = argparse.ArgumentParser("""QUAST - quality control""")
        parser.add_argument('--quast'        , help="compute assembly stats using quast", required=False, action="store_true")
        parser.add_argument('--contigs', nargs='*', help="input contigs (multiple assemblies are allowed)", required=True)
        parser.add_argument('--ref',required = False,type=str, help="reference fasta")
        parser.add_argument('--output',required = True,type=str, help="output folder")
	parser.add_argument('--features',nargs='*',required = False,type=str, help="Feature BED/GFF file")
        parser.add_argument('--len',default=100,type=int, help="minimum contig length (default= 100 bp)")
        args= parser.parse_args()

	quast="quast.py {} --output-dir {} --min-contig {}".format(" ".join(args.contigs),args.output,args.len)

	if args.ref:
		quast+=" -r {}".format(args.ref)

	if args.features:
		quast+=" -g {}".format(" ".join(args.features))
	os.system(quast)

elif args.align:
	parser = argparse.ArgumentParser("""Assemblatron align - align contigs to the reference using bwa mem""")
	parser.add_argument('--align'          , help="align contigs to reference using bwa mem", required=False, action="store_true")
	parser.add_argument('--ref',required = True,type=str, help="reference fasta")
	parser.add_argument('--mem'      , help="maximum  mempory per thread (gigabytes)", type=int, default=2)
	parser.add_argument('--cores'       ,type=int, default = 8, help="number of cores (default = 2)", required=False)
	parser.add_argument('--contigs',required = True,type=str, help="input contigs")
	parser.add_argument('--prefix',required = True,type=str, help="output prefix")
	args= parser.parse_args()

	os.system("bwa mem -x intractg -t {} {} {} | samtools view -Sbh - | samtools sort -m {}G - > {}.bam".format(args.cores,args.ref,args.contigs,args.mem,args.prefix))
	os.system( "samtools index {}.bam".format(args.prefix) )
	
elif args.fasta:

	parser = argparse.ArgumentParser("""Assemblatron fasta - converts bam to fasta using samtools""")
	parser.add_argument('--fasta'      , help="convert bam to fasta", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	args= parser.parse_args()

	os.system("samtools view -h -F 2048 {} | samtools view -F 1024 -Sh - | samtools view -Subh -F 256 - | samtools fasta - ".format(args.bam))

elif args.scaffold:
	parser = argparse.ArgumentParser("""Assemblatron scaffold - Perform scaffolding using BESST""")
	parser.add_argument('--scaffold'      , help="perform scaffolding using BESST", required=False, action="store_true")
	parser.add_argument('--fastq',required = True,type=str, help="input paired reads")
	parser.add_argument('--contigs',required = True,type=str, help="input contigs (fasta format)")
	parser.add_argument('--output',required = True,type=str, help="output folder")
	parser.add_argument('--rf'      , help="reverse forward orientation (i.e mate pair)", required=False, action="store_true")
	parser.add_argument('--fr'      , help="forward reverse orientation (i.e standard paired reads, default setting)", required=False, action="store_true")
	parser.add_argument('--tmpdir'      , help="write reads-to-contig bam to $TMPDIR ", required=False, action="store_true")
	parser.add_argument('--filename',required = True, help="filename of the output files (default = same as the contigs)")
	parser.add_argument('--mem'      , help="maximum  mempry per thread (gigabytes)", type=int, default=4)
        parser.add_argument('--iter'      , help="Number of itterations (default = 500000)", type=int, default=500000)
	parser.add_argument('--cores'       ,type=int, default = 8, help="number of cores (default = 2)", required=False)
        parser.add_argument('-q'       ,type=int, help="minimum mapping quality for scaffolding", required=False)
        parser.add_argument('-p'       ,type=int, help="minimum number of read-pairs to create edge", required=False)
	args= parser.parse_args()

	args.prefix=args.filename

	if not args.prefix:
		args.prefix=contigs.split("/")[-1].split(".")[0]
	if args.tmpdir:
		args.bam="$TMPDIR/{}.bam".format(args.prefix)
	else:
		args.bam="{}/{}.bam".format(args.output,args.prefix)

	os.system("mkdir -p {}".format(args.output))
	os.system("bwa index {}".format(args.contigs))
	os.system("bwa mem -p -t {} {} {} | samtools view -Sbh - | samtools sort -@ {} -m {}G - > {}".format(args.cores,args.contigs,args.fastq,args.cores,args.mem,args.bam))
	os.system("samtools index {}".format(args.bam))

	if args.rf:
		besst="runBESST -c {} -f {} -orientation rf -o {} -plots --iter {}".format(args.contigs,args.bam,args.output,args.iter)
	else:
		besst="runBESST -c {} -f {} -orientation fr -o {} -plots --iter {}".format(args.contigs,args.bam,args.output,args.iter)

	if args.q:
		besst+= " --min_mapq {}".format(args.q)

	if args.p:
                besst+= " -e {}".format(args.p)

	os.system(besst)
elif args.fastq:

	parser = argparse.ArgumentParser("""Assemblatron fastq - converts bam to fastq using samtools""")
	parser.add_argument('--fastq'      , help="convert bam to fastq", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (aligned reads)")
	parser.add_argument('--cores'       ,type=int, default = 8, help="number of cores (default = 2)", required=False)
	parser.add_argument('--sort'      , help="sort the  output fastq based on the read names(required for later aligning paired-reads)", required=False, action="store_true")
	parser.add_argument('--mem'      , help="maximum  mempry per thread (gigabytes)", type=int, default=4)
	args= parser.parse_args()

	if args.sort:
		os.system("samtools view -h -F 2048 {} | samtools view -F 1024 -Sh - | samtools view -Subh -F 256 - | samtools sort -n -m {}G -@ {} - | samtools bam2fq -".format(args.bam,args.mem,args.cores))
	else:
		os.system("samtools view -h -F 2048 {} | samtools view -F 1024 -Sh - | samtools view -Subh -F 256 - | samtools bam2fq -".format(args.bam))

elif args.snv:
	parser = argparse.ArgumentParser("""Assemblatron snv - call snvs using samtools pileup and bcftools""")
	parser.add_argument('--snv'             , help="call snv from the aligned contigs", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	parser.add_argument('--ref',required = True,type=str, help="reference fasta")

	args= parser.parse_args()
	os.system("{}/fermikit/htsbox/htsbox pileup -cuf {} {} ".format(wd,args.ref,args.bam))
else:
	parser.print_help()
