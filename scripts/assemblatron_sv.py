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

def report_splits(bam,working_dir,args,wd):
    bam_prefix=args.bam.split("/")[-1]
    prefix=args.working_dir + "/" +  bam_prefix[0:-4]
    os.system("{}/htsbox/htsbox abreak -b {}_splits.bam -c -f {} > {}_raw.vcf".format(wd,prefix,args.ref,prefix))
    contig_length={}
    mapQ={}

    cigarOP={}
    cigarList=[]
    with os.popen("samtools view {}_splits.bam".format(prefix)) as pipe:
        for line in pipe:
            content=line.strip().split()
            if not content[0] in contig_length:
                contig_length[content[0]]=0
            length=len(content[9])
            if length > contig_length[content[0]]:
                contig_length[content[0]]=length
            if not content[0] in cigarOP:
                cigarOP[content[0]] = []
            else:
                cigarList.append(content[5])
    

    output=[]
    i=0
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
        content=line.strip().split()
        info=content[7]+";CIGAR={}".format(cigarList[i])
        content[7] = info
        output.append("\t".join(content) + "\n")        
        i+=1

    coverage_data=coverage(args)
    f=open("{}_raw.vcf".format(prefix),"w")
    for line in output:
        if line[0] == "#":
            if "##FILTER" in line:
                f.write("##FILTER=<ID=MinQ,Description=\"minimum quality less than {}\">\n".format(args.q))
                f.write("##FILTER=<ID=HighCoverage,Description=\"coverage exceeding {}\">\n".format(args.max_coverage))
            if "##INFO=<ID=SVTYPE," in line:
                f.write("##INFO=<ID=CIGAR,Number=1,Type=String,Description=\"The cigar operation\">\n")
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
    os.system("samtools view -q {} -h -F 4 {} | grep -E \"@|SA:\"  | samtools view -Shb -@ {} - | samtools sort -n -@ {} - > {}_splits.bam".format(args.q,bam,args.cores,args.cores,prefix))
    try:
        if os.path.getsize( '{}_splits.bam'.format(prefix) ):
           found = True
    except:
        found= False 
    return found
    
def compute_coverage(args,wd):
    bam_prefix=args.bam.split("/")[-1]
    prefix=args.working_dir + "/" +  bam_prefix[0:-4]
    os.system("{}/TIDDIT/bin/TIDDIT --cov -b {} -o {} -z 100 ".format(wd,args.bam,prefix))


def main(args,assemblatron_home):
    if not os.path.isdir( args.working_dir ):
        os.makedirs( args.working_dir )
    found=find_splits(args.bam,args.working_dir,args)
    #found=True
    if found:
        compute_coverage(args,assemblatron_home)
        report_splits(args.bam,args.working_dir,args,assemblatron_home)
    else:
        print "warning: found no supplementary alignments!"
