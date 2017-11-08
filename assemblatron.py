import os
import sys
import numpy
import argparse
#import readVCF
import itertools

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

    f=open("{}_raw.vcf".format(prefix),"w")
    for line in output:
        f.write(line)

    os.system("svdb --merge --vcf {}_raw.vcf --bnd_distance 1 --overlap 1 --same_order > {}_merged.vcf".format(prefix,prefix) )
    os.system("vcf-sort -c {}_merged.vcf".format(prefix))

def find_splits(bam,working_dir,args):
    bam_prefix=bam.split("/")[-1]
    prefix=working_dir + "/" +  bam_prefix[0:-4]
    found = False
    os.system("samtools view -h -q {} -F 4 {} | grep -E \"@|SA:\"  > {}_splits.sam".format(args.q,bam,prefix))
    output=[]
    for line in open('{}_splits.sam'.format(prefix)):
        if line[0] == "@":
            output.append(line)
            continue
        content=line.strip().split()
        if len(content[9]) > args.len_aln:
            output.append(line)


    f=open('{}_splits.sam'.format(prefix),"w")
    for line in output:
        f.write(line)
    f.close()

    os.system("samtools view -h -Shb {}_splits.sam | samtools sort -n - {}_splits".format(prefix,prefix))
    if os.path.getsize( '{}_splits.sam'.format(prefix) ):
        found = True

    return found
    

def main(args):
    if not os.path.isdir( args.working_dir ):
        os.makedirs( args.working_dir )
    #found=find_splits(args.bam,args.working_dir,args)
    found=True
    if found:
        #clean_splits(args)
        report_splits(args.bam,args.working_dir,args)
    else:
        print "warning: found no supplementary alignments!"




parser = argparse.ArgumentParser("""Assemblatron - a variant caller using aligned contigs""")
parser.add_argument('--working_dir',default="work",type=str, help="temporary analysis files will be stored here(default=work)")
parser.add_argument('--bam',required = True,type=str, help="input bam")
parser.add_argument('--ref',required = True,type=str, help="reference fasta")
parser.add_argument('--q',type=int, default =0, help="minimum allowed mapping quality(default = 0)", required=False)
parser.add_argument('--len_aln'       ,type=int, default = 20, help="minimum length of alignments(default = 20)", required=False)
parser.add_argument('--len_ctg'       ,type=int, default = 1000, help="minimum contig length(default = 1000)", required=False)
parser.add_argument('--min_size'       ,type=int, default = 100, help="minimum variant size)", required=False)  
args= parser.parse_args()
main(args)

