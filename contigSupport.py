import os
import sys
import argparse
import numpy

wd= os.path.dirname(os.path.realpath(__file__))
def main(args,wd):
    os.system("mkdir {}".format(args.wd))
    os.system("samtools bam2fq {} | {}/seqtk/seqtk seq -a > {}".format(args.contigs,wd,args.contigs.replace(".bam",".fa")))
    args.contigs=args.contigs.replace(".bam",".fa")
    os.system("bwa index {}".format(args.contigs))

    #compute the coverage across the input bam and the contigs
    
    #print "samtools bam2fq {} | bwa mem -t {} {} -p - | samtools view -Sbh - | samtools sort - {}/{}".format(args.bam,args.cores,args.contigs,args.wd,args.prefix)
    os.system("samtools sort -t {} -n {} {}/{}.sorted".format(args.cores,args.bam,args.wd,args.prefix))
    os.system("samtools bam2fq {}/{}.sorted.bam | bwa mem -t {} {} -p - | samtools view -Sbh - | samtools sort - {}/{}".format(args.wd,args.prefix,args.cores,args.contigs,args.wd,args.prefix))
    os.system("{}/TIDDIT/bin/TIDDIT --cov -b {} -z 10000 -o {}/{}".format(wd,args.bam,args.wd,args.prefix+".mapped_coverage"))
    os.system("{}/TIDDIT/bin/TIDDIT --cov -b {}/{}.bam -z 100 -o {}/{}".format(wd,args.wd,args.prefix,args.wd,args.prefix+".contig_coverage"))

    genomic_coverage=[]
    for line in open("{}/{}.mapped_coverage.tab".format(args.wd,args.prefix)):
        if line[0] == "#":
            continue

        content=line.strip().split()
        if float(content[3]) > 0:
            genomic_coverage.append( float(content[3]) )        
        
    avg_coverage=numpy.average(genomic_coverage)
    del genomic_coverage

    reads=10000000
    limit=100000
    i=0
    insert_sizes=[]
    with os.popen("samtools view {}".format(args.bam)) as pipe:
        for line in pipe:
            content=line.strip().split()
            ins=abs(int(content[8]))
            if ins < limit and content[6] == "=":
                insert_sizes.append(ins)
            i+=1
            if i > reads:
                break

    mean_ins=numpy.average(insert_sizes)
    std_ins=numpy.std(insert_sizes)
    del insert_sizes

    contig_stats={}
    #compute the number of read pairs split between two contigs, and the number of read-pairs mapping to the same contig
    unmapped=0
    mapped=0
    links={}
    with os.popen("samtools view {}/{}.bam".format(args.wd,args.prefix)) as pipe:
        for line in pipe:
            content=line.strip().split()
            if content[6] == "*":
                unmapped+=1
                continue
            mapped+=1

            if not content[2] in contig_stats:
                contig_stats[content[2]] = [0,0,0]
                links[content[2]]=set([])
            if content[6] == "=":
                if abs( abs(int(content[8]))-mean_ins) <= std_ins*3:
                    contig_stats[content[2]][1] +=1
                else:
                    contig_stats[content[2]][0] +=1
            else:
                contig_stats[content[2]][2] +=1
                links[content[2]].add(content[6])


    #compute the average content across each contig, as well as the average coverage across the regions were the contigs map
    contig_coverage={}
    for line in open("{}/{}.contig_coverage.tab".format(args.wd,args.prefix)):
        if line[0] == "#":
            continue
        content=line.strip().split()
        if not content[0] in contig_coverage:
            contig_coverage[content[0]]=[]
        contig_coverage[content[0]].append(float(content[3])) 

    for contig in contig_coverage:
        contig_coverage[contig]=numpy.average(contig_coverage[contig])

    contig_len={}
    with os.popen("samtools view -H {}/{}.bam".format(args.wd,args.prefix)) as pipe:
        for line in pipe:
            if "@SQ" in line:
                content=line.strip().split()
                contig_len[content[1].split("SN:")[-1]]=content[2].split("LN:")[-1]


    #print summary
    print ("#library_stats\taverage_coverage={},mean_ins={},std_ins={},unmapped={},mapped={}".format(avg_coverage,mean_ins,std_ins,unmapped,mapped))
    print ("#contig_id\tlength\tcontig_coverage\tdiscordant_pairs\tconcordant_pairs\tinter_contig_pairs\ttotal_pairs\tlinks_to_contigs")
    for contig in contig_stats:
        if contig == "*":
            continue
        print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(contig,contig_len[contig],contig_coverage[contig],contig_stats[contig][0],contig_stats[contig][1],contig_stats[contig][2],sum(contig_stats[contig]),",".join(links[contig]) )

    return()

parser = argparse.ArgumentParser("""contigSupport - quality control of contigs""")
parser.add_argument('--contigs',type=str,required=True, help="aligned contigs produced through de novo assembly")
parser.add_argument('--bam',type=str,required=True, help="bam file containing the reads")
parser.add_argument('--prefix',required = True,type=str, help="prefix of the output files")
parser.add_argument('--wd',type=str,default="contigAnalysis",help="the files will be  printed to this directory")
parser.add_argument('--cores',type=int, default =16, help="number of cores (default = 16)")
args= parser.parse_args()
main(args,wd)
