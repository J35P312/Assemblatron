import os
import sys
import math
import numpy
import argparse
#import readVCF
import itertools

def coverage(args):
    coverage_data={}
    bam_prefix=args.bam.split("/")[-1]
    prefix=args.working_dir + "/" +  bam_prefix[0:-4]
    for line in open(prefix+".tab"):
        if line[0] == "#":
            continue
        content=line.strip().split()
        if not content[0] in coverage_data:
            coverage_data[content[0]]=[]
        coverage_data[content[0]].append(float(content[3]))

    for chromosome in coverage_data:
        coverage_data[chromosome]=numpy.array(coverage_data[chromosome])

    return(coverage_data)

def report_splits(bam,working_dir,args):
    bam_prefix=args.bam.split("/")[-1]
    prefix=args.working_dir + "/" +  bam_prefix[0:-4]
    wd= os.path.dirname(os.path.realpath(__file__))
    os.system("{}/htsbox/htsbox abreak -b {}_splits.bam -c -f {} > {}_raw.vcf".format(wd,prefix,args.ref,prefix))
    contig_length={}
    mapQ={}
    for line in open('{}_splits.sam'.format(prefix)):
        if line[0] == "@":
            continue
        content=line.strip().split()
        if not content[0] in contig_length:
            contig_length[content[0]]=0
        length=len(content[9])
        if length > contig_length[content[0]]:
            contig_length[content[0]]=length

    

    output=[]
    for line in open("{}_raw.vcf".format(prefix)):
        if line[0] == "#":
            output.append(line)
            continue
        contigid=line.split("CONTIGID=")[-1]
        if contig_length[contigid.strip()] < args.len_ctg:
            continue
        if "SVLEN=" in line:
            if args.min_size > int(line.split("SVLEN=")[-1].split(";")[0]):
                continue
        output.append(line)        

    coverage_data=coverage(args)
    f=open("{}_raw.vcf".format(prefix),"w")
    for line in output:
        if line[0] == "#":
            if "##FILTER" in line:
                f.write("##FILTER=<ID=MinQ,Description=\"minimum quality less than {}\">\n".format(args.q))
                f.write("##FILTER=<ID=HighCoverage,Description=\"coverage exceeding {}\">\n".format(args.max_coverage))
            else:
                f.write(line)
            continue
        content=line.strip().split()
        chrA=content[0]
        posA=int(math.floor(int(content[1])/100.0))

        rA= coverage_data[chrA][posA-1:posA+2]
        quality=int(line.split("MINMAPQ=")[-1].split(";")[0])
        if len(rA):
            if max(rA) > args.max_coverage:
                content[6]="HighCoverage"
            elif quality < args.q:
                content[6]="MinQ"                
            else:
                content[6]="PASS"          

        if "<" in content[4]:
            chrB=chrA
            posB=int(math.floor(int( content[7].split(";END=")[-1].split(";")[0] )/100.0))
        elif ":" in content[4]:
            chrB=content[4].split(":")[0].split("[")[-1].split("]")[-1]
            posB=int(math.floor(int(content[4].split(":")[-1].split("[")[0].split("]")[0])/100.0))

        rB=coverage_data[chrB][posB-1:posB+2]
        if len(rB):
            if max(rB) > args.max_coverage:
                content[6]="HighCoverage"

        f.write("\t".join(content)+"\n")
    f.close()

    os.system("svdb --merge --vcf {}_raw.vcf --bnd_distance 1 --overlap 1 --same_order > {}_merged.vcf".format(prefix,prefix) )
    os.system("vcf-sort -c {}_merged.vcf".format(prefix))    

def find_splits(bam,working_dir,args):
    bam_prefix=bam.split("/")[-1]
    prefix=working_dir + "/" +  bam_prefix[0:-4]
    found = False
    os.system("samtools view -h -F 4 {} | grep -E \"@|SA:\"  > {}_splits.sam".format(args.q,bam,prefix))
    output=[]
    for line in open('{}_splits.sam'.format(prefix)):
        if line[0] == "@":
            output.append(line)
            continue
        content=line.strip().split()
        output.append(line)


    f=open('{}_splits.sam'.format(prefix),"w")
    for line in output:
        f.write(line)
    f.close()

    os.system("samtools view -h -Shb {}_splits.sam | samtools sort -n - {}_splits".format(prefix,prefix))
    if os.path.getsize( '{}_splits.sam'.format(prefix) ):
        found = True

    return found
    
def compute_coverage(args):
    wd= os.path.dirname(os.path.realpath(__file__))
    bam_prefix=args.bam.split("/")[-1]
    prefix=args.working_dir + "/" +  bam_prefix[0:-4]
    #os.system("{}/TIDDIT/bin/TIDDIT --cov -b {} -o {} -z 100 ".format(wd,args.bam,prefix))


def main(args):
    if not os.path.isdir( args.working_dir ):
        os.makedirs( args.working_dir )
    #found=find_splits(args.bam,args.working_dir,args)
    found=True
    if found:
        compute_coverage(args)
        report_splits(args.bam,args.working_dir,args)
    else:
        print "warning: found no supplementary alignments!"




parser = argparse.ArgumentParser("""Assemblatron - a variant caller using aligned contigs""")
parser.add_argument('--working_dir',default="work",type=str, help="temporary analysis files will be stored here(default=work)")
parser.add_argument('--bam',required = True,type=str, help="input bam")
parser.add_argument('--ref',required = True,type=str, help="reference fasta")
parser.add_argument('--q',type=int, default =10 ,help="minimum allowed mapping quality(default = 10)", required=False)
parser.add_argument('--len_ctg'       ,type=int, default = 1000, help="minimum contig length(default = 1000)", required=False)
parser.add_argument('--max_coverage'       ,type=int, default = 8, help="calls from regions exceeding the maximum coverage are filtered", required=False)
parser.add_argument('--min_size'       ,type=int, default = 100, help="minimum variant size)", required=False)  
args= parser.parse_args()
main(args)

