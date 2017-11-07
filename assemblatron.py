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
    os.system("{}/htsbox/htsbox abreak {}_splits.sam -c".format(wd,prefix))

def find_splits(bam,working_dir,args):
    bam_prefix=bam.split("/")[-1]
    prefix=working_dir + "/" +  bam_prefix[0:-4]
    found = False
    #print "samtools view {} | grep \"SA:\" | grep -E :{},|;{}, > test.sam".format(bam)
    os.system("samtools view -h -q {} -F 4 {} | grep -E \"@|SA:\"  > {}_splits.sam".format(args.q,bam,prefix))
    if os.path.getsize( '{}_splits.sam'.format(prefix) ):
        found = True

    return found

def clean_splits(args):
    bam_prefix=args.bam.split("/")[-1]
    prefix=args.working_dir + "/" +  bam_prefix[0:-4]
    first=True
    lines=[]
    for line in open("{}_splits.sam".format(prefix)):
        content=line.strip().split()
        if "@" in line:
            lines.append(line)
            continue
        if first:
            line_id=int(content[0])
            sequence=content[9]
            lines.append(line)
            first = False
            continue
        if int(content[0]) == line_id+1 or int(content[0]) == line_id-1 and sequence == content[9]:
            continue
        if len(sequence) < args.len_ctg:
            continue
        lines.append(line)
        line_id=int(content[0])
        sequence=content[0]
    f=open("{}_splits.sam".format(prefix),"w")
    for line in lines:
        f.write(line)
    f.close()


def main(args):
    if not os.path.isdir( args.working_dir ):
        os.makedirs( args.working_dir )
    found=find_splits(args.bam,args.working_dir,args)
    found = True
    if found:
        report_splits(args.bam,args.working_dir,args)
    else:
        print "warning: found no supplementary alignments!"




parser = argparse.ArgumentParser("""Assemblatron - a variant caller using aligned contigs""")
parser.add_argument('--working_dir',default="work",type=str, help="temporary analysis files will be stored here(default=work)")
parser.add_argument('--bam',required = True,type=str, help="input bam")
parser.add_argument('--q',type=int, default =10, help="minimum allowed mapping quality(default = 10)", required=False)
parser.add_argument('--len'       ,type=int, default = 50, help="minimum length of alignments(default = 50)", required=False)
parser.add_argument('--len_ctg'       ,type=int, default = 1000, help="minimum contig length(default = 1000)", required=False) 
args= parser.parse_args()
main(args)

