import numpy
import math
import itertools
import os

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
	unmapped_contigs=0
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
					c=numpy.zeros(bins)
					chromosome_order.append(chromosome)

	i=0
	with os.popen("samtools view {}".format(args.bam)) as pipe:
		for line in pipe:
			i+=1
			content=line.strip().split()
			if content[2] != "*":

				sequence_length=compute_aln_length(content[5])
				binSize=50
				element=int(math.floor(int(content[3])/50.0))
				q=int(content[4])
				chromosome=content[2]
                #coverage_structure[content[2]][pos]

				filled_elements=int(math.floor(sequence_length/50.0))
				remainingRead=(sequence_length % binSize)/50
				for i in range(0,filled_elements):
					coverage_structure[content[2]][element]+=1;
					if q > 30:
						coverage_structure_20[content[2]][element]+=1;
					element+=1

				if (remainingRead > 0 and len(coverage_structure[content[2]]) > element+1):
					coverage_structure[content[2]][element]+=remainingRead;
					if q > 30:
						coverage_structure_20[content[2]][element]+=remainingRead;    

			else:
				unmapped+=len(content[9])
				unmapped_contigs+=1

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

	for contig in sorted(contig_sizes, reverse=True):
		n+= contig
	
		if n <= assembly_size/2: 
			N_50= contig
			L_50+=1

		if n <= assembly_size*0.9:  
			N_90= contig
			L_90+=1      
			
            

	print "\t".join	(["file","N50","L50","N90","L90","number_of_contigs","assembly_size","longest_contig","zero_coverage","covered_bases","zero_coverage(q>30)","covered_bases(q>30)","reference_length","unmapped_bases","unmapped%","unmapped_contigs"])
	print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(args.bam,N_50,L_50,N_90,L_90,len(contig_sizes),assembly_size,max(contig_sizes),uncovered/float(total),covered/float(total),uncovered_20/float(total),covered_20/float(total),total,unmapped,unmapped/float(assembly_size),unmapped_contigs)
