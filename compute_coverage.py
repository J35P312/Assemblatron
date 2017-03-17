import sys
import numpy
#first argument input coverage file from TIDDIT --cov
bin_size=50
genome=[]
for line in open(sys.argv[1]):
    if line[0] == "#":
        continue
    if 50 > abs(int(line.strip().split("\t")[1]) - int(line.strip().split("\t")[2])): 
        continue
    genome.append( float(line.strip().split("\t")[3]) )


genome=numpy.array(genome)

coverage_values=range(0,101)
total_bins=len(genome)

print "#percentage of genome zero coverage"
print float(len(genome[numpy.where(genome  ==  0)]))/float(total_bins)
print "#cumulative_coverage\tpercentage"
for val in coverage_values:
    print( "{}\t{}".format(val, float(len( genome[ numpy.where(genome >= val) ]))/total_bins ) )


