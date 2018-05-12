import sys
import os
import argparse
import numpy
import math
import itertools

wd=os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, '{}/scripts'.format(wd))
import call

def compute_aln_length(cigar):
    length=0
    #compute the length of the contig
    SC = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
    for i in range(0,len(SC)/2):
        if SC[i*2+1] == "M":
            length += int( SC[i*2] )
    return length

def assembly_stats(args):
	uncovered=0
	covered=0
	uncovered_20=0
	covered_20=0
	total=0
	unmapped=0

	coverage_structure={}
	coverage_structure_20={}
	chromosome_order=[]
	with os.popen("samtools view -H {}".format(args.bam)) as pipe:
		for line in pipe:
			if line[0] == "@":
				if "SN:" in line:
					content=line.strip().split()
					chromosome=content[1].split("SN:")[-1]
					length=int(content[2].split("LN:")[-1])
					bins=int( math.ceil(length/float(50)) )
					coverage_structure[chromosome]=numpy.zeros(bins)
					coverage_structure_20[chromosome]=numpy.zeros(bins)
					chromosome_order.append(chromosome)

	with os.popen("samtools view {}".format(args.bam)) as pipe:
		for line in pipe:
			content=line.strip().split()
			if content[2] != "*":
				q=int(content[4])
				aln_len=compute_aln_length(content[5])

				if content[2] in chromosome_order:
					pos=int(math.floor(int(content[3])/50))
					to_add=int(math.floor(aln_len/50.0))
					for i in range(0,to_add):
						if i < to_add-1:
							coverage_structure[content[2]][pos]+= 1
							if q > 30:
								coverage_structure_20[content[2]][pos]+= 1
						else:
							coverage_structure[content[2]][pos]+=(len(content[9])-(to_add-1)*50)/50.0
							if q > 30:
								coverage_structure[content[2]][pos]+=(len(content[9])-(to_add-1)*50)/50.0
						pos+=1

			else:
				unmapped+=len(content[9])

	for chromosome in coverage_structure:
		for i in range(0,len(coverage_structure[chromosome])):
			if coverage_structure[chromosome][i] >= 1:
				covered += 50
			else:
				uncovered += 50

			if coverage_structure_20[chromosome][i] >= 1:
				covered_20 += 50
			else:
				uncovered_20 += 50

			total += 50
        

	contig_sizes=[]
	contig_size=0

	with os.popen("samtools view -h -F 2048 {} | samtools view -F 1024 -Sh - | samtools view -Subh -F 256 | samtools fasta -".format(args.bam)) as pipe:
		for line in pipe:
			if ">" == line[0]:
				if contig_size:
					contig_sizes.append(contig_size)
					contig_size=0
				continue
			contig_size += len(line.strip())

	assembly_size=sum(contig_sizes)
	n=0
	N_50=0
	L_50=0
	N_90=0
	L_90=0

	for contig in sorted(contig_sizes):
		n+= contig
	
		if n <= assembly_size/2: 
			N_50= contig
			L_50+=1

		if n <= assembly_size*0.9:  
			N_90= contig
			L_90+=1      
			
            

	print "\t".join	(["file","N50","L50","N90","L90","number_of_contigs","assembly_size","longest_contig","zero_coverage","covered_bases","zero_coverage(q>30)","covered_bases(q>30)","reference_length","unmapped_bases","unmapped%"])
	print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(args.bam,N_50,L_50,N_90,L_90,len(contig_sizes),assembly_size,max(contig_sizes),uncovered/float(total),covered/float(total),uncovered_20/float(total),covered_20/float(total),total,unmapped,unmapped/float(assembly_size))


def assemble(args,wd):
	fermi="{}/fermikit/fermi.kit/".format(wd)

	#apply bloom filter and build the index
	ropebwt="{}/ropebwt2 -m {} -dNCr - > {}.fmd 2> {}.fmd.log".format(fermi,args.batch,args.prefix,args.prefix)
	bfc="{}/bfc -1s {} -k {} -t {} {} 2> {}.flt.fq.gz.log".format(fermi,args.z,args.l,args.cores,args.fastq,args.prefix)
	os.system("{} | {}".format(bfc,ropebwt))
	#assemble
	os.system( "{}/fermi2 assemble -l {} -m {} -t {} {}.fmd 2> {}.pre.gz.log | gzip -1 > {}.pre.gz".format(fermi,args.l,args.m,args.cores,args.prefix,args.prefix,args.prefix) )
	os.system("{}/fermi2 simplify -CSo 66 -m {} -T 61 {}.pre.gz 2>  {}.mag.gz.log > {}.fasta".format(fermi,args.m,args.prefix,args.prefix, args.prefix))

	if args.align:
		os.system("bwa mem -x intractg -t {} {} {}.fasta | samtools view -Sbh - | samtools sort -m 2G - > {}.bam".format(args.cores,args.ref,args.prefix,args.prefix))
		os.system("samtools index {}.bam".format(args.prefix) )

version = "0.0.0"
parser = argparse.ArgumentParser("""Assemblatron: a de novo assembly  pipeline""".format(version),add_help=False)
parser.add_argument('--assemble'       , help="Perform de novo assembly using the Fermi2 assembler", required=False, action="store_true")
parser.add_argument('--scaffold'      , help="perform scaffolding using BESST", required=False, action="store_true")
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
	parser.add_argument('--cores',type=int, default =16, help="number of cores (default = 16)")
	parser.add_argument('--batch',type=str, default ="20g", help="batch size for multi-string indexing; 0 for single-string (default=20g)")
	parser.add_argument('-z',type=str, default ="3G", help="genome size (use K,M, or G) (default = 3G)")
	parser.add_argument('-l',type=int, default =81, help="min match (default = 81)")
	parser.add_argument('-m',type=int, default =100, help="min merge length (default = 100)")
	parser.add_argument('--align', help="align contigs to reference using bwa mem", required=False, action="store_true")
	parser.add_argument('--ref',type=str, help="reference fasta, required for alignment of the contigs")
	args= parser.parse_args()

	if not args.align or (args.align and args.ref):
		assemble(args,wd)
	elif args.align and not args.ref:
		print ("you need a reference to align the contigs: Please supply the reference path through the --ref parameter")

elif args.sv:

	parser = argparse.ArgumentParser("""Assemblatron sv - a variant caller using aligned contigs""")
	parser.add_argument('--sv'        , help="call SV from the aligned contigs", required=False, action="store_true")
	parser.add_argument('--bam',required = True,type=str, help="input bam (contigs)")
	parser.add_argument('--q',type=int, default =10 ,help="minimum allowed mapping quality(default = 10)", required=False)
	parser.add_argument('--len_ctg'       ,type=int, default = 35, help="minimum uniqyelly mapped contig length(default = 35)", required=False)
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

	assembly_stats(args)

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
		os.system("runBESST -c {} -f {} -orientation rf -o {} --iter {}".format(args.contigs,args.bam,args.output,args.iter))
	else:
		os.system("runBESST -c {} -f {} -orientation fr -o {} --iter {}".format(args.contigs,args.bam,args.output,args.iter))
	
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
	os.system("samtools mpileup -uf {} {}  | bcftools call -cv".format(args.ref,args.bam))
else:
	parser.print_help()
