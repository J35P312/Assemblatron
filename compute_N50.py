import sys
# the first argument is the input fa file, the software will compute the N50, number of contigs, longest contig, and total assembly size to stdout


contig_sizes=[]
contig_size=0
for line in open(sys.argv[1]):
    if ">" == line[0]:
       if contig_size:
          contig_sizes.append(contig_size)
          contig_size=0
       continue
    contig_size += len(line.strip())


assembly_size=sum(contig_sizes)
n=0
N=0
for contig in sorted(contig_sizes):
    n+= contig
    if n >= assembly_size/2: 
       N= contig
       break

print "\t".join(["file","N50","number_of_contigs","assembly_size","longest_contig"])
print "{}\t{}\t{}\t{}\t{}".format(sys.argv[1],N,len(contig_sizes),assembly_size,max(contig_sizes))

