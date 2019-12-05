#!/usr/bin/env python
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
	if not args.noec:
		#run kmc
		kmc_prefix="{}.kmc".format(args.prefix)
		kmc="kmc -k{} -m30 {} {} {}".format(args.k,args.fastq,kmc_prefix,args.tmp)
		os.system(kmc)

		#apply bloom filter and build the index
		ropebwt="{}/ropebwt2 -m {} -dNCr - > {}.fmd 2> {}.fmd.log".format(fermi,args.batch,args.prefix,args.prefix)
		bfc="{}/bfc-kmc -k {} -T -t {} {} {} 2> {}.flt.fq.gz.log".format(bfc_path,args.k,args.cores,kmc_prefix,args.fastq,args.prefix)
		os.system("{} | {}".format(bfc,ropebwt))
	else:
		ropebwt="{}/ropebwt2 -m {} -dNCr {} > {}.fmd".format(fermi,args.batch,args.fastq,args.prefix)
		os.system(ropebwt)

	#assemble
	os.system( "{}/fermi2 assemble -l {} -t {} {}.fmd | gzip -1 > {}.pre.gz".format(fermi,args.l,args.cores,args.prefix,args.prefix,args.prefix) )
	if args.A:
		os.system("{}/fermi2 simplify -CA -R {} -d {} {}.pre.gz 2>  {}.mag.gz.log > {}.mag".format(fermi,args.r,args.r,args.prefix,args.prefix, args.prefix))
	else:
		os.system("{}/fermi2 simplify -CS -R {} -d {} {}.pre.gz 2>  {}.mag.gz.log > {}.mag".format(fermi,args.r,args.r,args.prefix,args.prefix, args.prefix))

	if args.align:
		os.system("bwa mem -x intractg -t {} {} {}.mag | samtools view -Sbh - | sambamba sort -m 10G -t /dev/stdin -o {}.bam".format(args.cores,args.ref,args.prefix,args.cores,args.prefix))

version = "1.0.0"
parser = argparse.ArgumentParser("""Assemblatron: a de novo assembly  pipeline""".format(version),add_help=False)
parser.add_argument('--assemble'       , help="Perform de novo assembly using the Fermi2 assembler", required=False, action="store_true")
parser.add_argument('--local'       , help="Perform local de novo assembly using the Fermi2 assembler", required=False, action="store_true")
parser.add_argument('--sv'             , help="call SV from the aligned contigs", required=False, action="store_true")
parser.add_argument('--snv'             , help="call snv from the aligned contigs", required=False, action="store_true")
parser.add_argument('--stats'          , help="compute assembly stats from aligned contigs", required=False, action="store_true")
parser.add_argument('--align'          , help="align contigs to reference using bwa mem", required=False, action="store_true")
parser.add_argument('--fasta'      , help="convert aligned contigs bam file to fasta", required=False, action="store_true")
parser.add_argument('--fastq'      , help="convert bam to fastq", required=False, action="store_true")

args, unknown = parser.parse_known_args()

if args.assemble:

	parser = argparse.ArgumentParser("""Assemblatron assemble - a wrapper for the fermi assembler""")
	parser.add_argument('--assemble'       , help="Perform de novo assembly using the Fermi2 assembler", required=False, action="store_true")
	parser.add_argument('--fastq',required = True, type=str, help="input fastq")
	parser.add_argument('--prefix',required = True,type=str, help="prefix of the output files")
	parser.add_argument('--noec'      , help="skip error correction, useful for small datasets (local assembly etc)", required=False, action="store_true")
	parser.add_argument('--cores',type=int, default =16, help="number of cores (default = 16)")
	parser.add_argument('--batch',type=str, default ="20g", help="batch size for multi-string indexing; 0 for single-string (default=20g)")
	parser.add_argument('-l',type=int, default =81, help="min match (default = 81)")
	parser.add_argument('-k',type=int, default =41, help="minimum kmer length for kmc/bfc error correction (default = 41)")
	parser.add_argument('-r',type=float, default =0.9, help="minimum coverlap ratio between vertices (default=0.9)")
	parser.add_argument('-A', help="Agressive bubble poping, produce more contigous assemblies, sometimes at cost of lower precisison",required=False, action="store_true")
	parser.add_argument('--align', help="align contigs to reference using bwa mem", required=False, action="store_true")
	parser.add_argument('--ref',type=str, help="reference fasta, required for alignment of the contigs")
	parser.add_argument('--tmp',type=str,default="$TMPDIR", help="tmp directory, kmc will write tmp files here (default=$TMPDIR)")
	args= parser.parse_args()

	if not os.path.isdir(args.tmp) and not args.noec:
		print "error: no such directory {}".format(args.tmp)
		print "set the --tmp variable to an existing folder"
		quit()

	if not args.align or (args.align and args.ref):
		assemble(args,wd)
	elif args.align and not args.ref:
		print ("you need a reference to align the contigs: Please supply the reference path through the --ref parameter")

elif args.local:
	parser = argparse.ArgumentParser("""Assemblatron local - a wrapper for the fermi assembler""")
	parser.add_argument('--local'       , help="Perform local de novo assembly using the Fermi2 assembler", required=False, action="store_true")
	parser.add_argument('--bam',required = True, type=str, help="input bam")
	parser.add_argument('-l',type=int, default =81, help="min match (default = 81)")
	parser.add_argument('--region',type=str,required=True, help="genomic region to assemble; using sthe ame format  as samtools. Multiple regions are allowed",nargs='*')
	parser.add_argument('-r',type=float, default =0.9, help="minimum coverlap ratio between vertices (default=0.9)")
	parser.add_argument('--ref',type=str, help="reference fasta, required for alignment of the contigs")
	args= parser.parse_args()


	fermi="{}/fermikit/fermi.kit/".format(wd)
	samtools="samtools view -bh {} {} | samtools fastq - ".format(args.bam," ".join(args.region)) 
	ropebwt="{}/ropebwt2 -dNCr - ".format(fermi)
	fermi_assemble="{}/fermi2 assemble -l {} - ".format(fermi,args.l) 
	fermi_simplify="{}/fermi2 simplify -C -R {} -d {} - ".format(fermi,args.r,args.r,args.l)

	os.system("{} | {} | {} | {}".format(samtools,ropebwt,fermi_assemble,fermi_simplify))

elif args.sv:

	parser = argparse.ArgumentParser("""Assemblatron sv - a variant caller using aligned contigs""")
	parser.add_argument('--sv'        , help="call SV from the aligned contigs", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	parser.add_argument('--q',type=int, default =10 ,help="minimum allowed mapping quality(default = 10)", required=False)
	parser.add_argument('--len_ctg'       ,type=int, default = 40, help="minimum uniquelly mapped contig length(default = 40)", required=False)
	parser.add_argument('--max_size'       ,type=int, help="filter variants exceeding the following value (bp) (default=infinity)", required=False)
	parser.add_argument('--max_contigs'       ,type=int, default = 8, help="filter breakpoint regions containing too many contings (default=8 contigs)", required=False)
	parser.add_argument('--max_coverage'       ,type=int, default = 3, help="filter breakpoint regions located in high coverage regions (default=3*avg_chromosomal_coverage)", required=False)
	parser.add_argument('--min_size'       ,type=int, default = 50, help="minimum variant size)", required=False)  
	parser.add_argument('--sample'       ,type=str, help="sample id, as shown in the format column", required=False)
	args= parser.parse_args()

	call.main(args)

elif args.stats:
	parser = argparse.ArgumentParser("""Assemblatron stats - compute assembly stats""")
	parser.add_argument('--stats'        , help="compute assembly stats from aligned contigs", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	args= parser.parse_args()

	stats.assembly_stats(args)

elif args.align:
	parser = argparse.ArgumentParser("""Assemblatron align - align contigs to the reference using bwa mem""")
	parser.add_argument('--align'          , help="align contigs to reference using bwa mem", required=False, action="store_true")
	parser.add_argument('--ref',required = True,type=str, help="reference fasta")
	parser.add_argument('--mem'      , help="maximum  memory (gigabytes)", type=int, default=10)
	parser.add_argument('--cores'       ,type=int, default = 8, help="number of cores (default = 8)", required=False)
	parser.add_argument('--contigs',required = True,type=str, help="input contigs")
	parser.add_argument('--prefix',required = True,type=str, help="output prefix")
	args= parser.parse_args()

	os.system("bwa mem -x intractg -t {} {} {} | samtools view -Sbh - | sambamba sort -m {}G /dev/stdin -o {}.bam".format(args.cores,args.ref,args.contigs,args.mem,args.prefix))
	os.system( "samtools index {}.bam".format(args.prefix) )
	
elif args.fasta:

	parser = argparse.ArgumentParser("""Assemblatron fasta - converts bam to fasta using samtools""")
	parser.add_argument('--fasta'      , help="convert bam to fasta", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	args= parser.parse_args()

	os.system("samtools fasta {} ".format(args.bam))

elif args.fastq:

	parser = argparse.ArgumentParser("""Assemblatron fastq - converts bam to fastq using samtools""")
	parser.add_argument('--fastq'      , help="convert bam to fastq", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (aligned reads)")
	parser.add_argument('--cores'       ,type=int, default = 8, help="number of cores (default = 2)", required=False)
	parser.add_argument('--sort'      , help="sort the  output fastq based on the read names(required for later aligning paired-reads)", required=False, action="store_true")
	parser.add_argument('--mem'      , help="maximum memory (gigabytes)", type=int, default=20)
	args= parser.parse_args()

	if args.sort:
		os.system("sambamba sort -n -m {}G -t {} {} -o /dev/stdout | samtools fastq -".format(args.mem,args.cores,args.bam))
	else:
		os.system("samtools fastq {}".format(args.bam))

elif args.snv:
	parser = argparse.ArgumentParser("""Assemblatron snv - call snvs and indels usisng htsbox pileup""")
	parser.add_argument('--snv'             , help="call snv from the aligned contigs", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	parser.add_argument('--ref',required = True,type=str, help="reference fasta")
	parser.add_argument('-q'       ,type=int,default=10, help="minimum mapping quality of contigs (default=10)")
	parser.add_argument('-r'       ,type=int,default=4, help="minimum number of reads supporting the SNV (default=4)")

	args= parser.parse_args()
	os.system("{}/fermikit/htsbox/htsbox pileup -d -c -V.05 -S 50 -q{} -s{} -f {} {} ".format(wd,args.q,args.r,args.ref,args.bam))
else:
	parser.print_help()
