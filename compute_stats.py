import sys
import os
# the first argument is the input bam file, the software will compute the N50, number of contigs, longest contig, and total assembly size to stdout

uncovered=0
covered=0
total=0
wd=os.path.dirname(os.path.realpath(__file__))

with os.popen("{}/TIDDIT/bin/TIDDIT --cov -b {} -z 100 -o stdout".format(wd,sys.argv[1])) as pipe:
    for line in pipe:
        if line[0] == "#":
            continue
        content= line.strip().split()
        if float(content[-2]) == 0:
           uncovered+=100
        else:
           covered+=100
        total+=100


contig_sizes=[]
contig_size=0

with os.popen("samtools view -h -F 2048 {} | samtools view -F 1024 -Sh - | samtools view -Subh -F 256 | samtools fasta -".format(sys.argv[1])) as pipe:
    for line in pipe:
       if ">" == line[0]:
          if contig_size:
               contig_sizes.append(contig_size)
               contig_size=0
          continue
       contig_size += len(line.strip())

assembly_size=sum(contig_sizes)
n=0
N=0
L=0
for contig in sorted(contig_sizes):
    n+= contig
    L+=1
    if n >= assembly_size/2: 
       N= contig
       break

print "\t".join(["file","N50","L50","number_of_contigs","assembly_size","longest_contig","zero_coverage","covered_bases","reference_length"])
print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sys.argv[1],N,L,len(contig_sizes),assembly_size,max(contig_sizes),uncovered/float(total),covered/float(total),total)

