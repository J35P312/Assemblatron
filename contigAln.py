import os
import sys
import argparse
import numpy

wd= os.path.dirname(os.path.realpath(__file__))

def main(args,wd):

    args.prefix=args.contigs.replace(".bam","").split("/")[-1]
    os.system("mkdir -p {}".format(args.wd))
    os.system("samtools fasta {} > {}/{}.fa".format(args.contigs,args.wd,args.prefix))
    os.system("bwa index {}/{}.fa".format(args.wd,args.prefix))
    args.contigs="{}/{}.fa".format(args.wd,args.prefix)

    #align the reads to the scaffold    
    os.system("samtools view -h -F 2048 {} | samtools view -F 1024 -Sh - | samtools view -Subh -F 256 - | samtools sort -m 4G -l 0 -@ {} -n - | samtools bam2fq - | bwa mem -t {} {} -p - | samtools view -Sbh - | samtools sort -m 10G - > {}/{}.bam".format(args.bam,args.cores,args.cores,args.contigs,args.wd,args.prefix))
    os.system("samtools index {}/{}.bam".format(args.wd,args.prefix))
    return()

parser = argparse.ArgumentParser("""contigSupport - quality control of contigs""")
parser.add_argument('--contigs',type=str,required=True, help="aligned contigs produced through de novo assembly")
parser.add_argument('--bam',type=str,required=True, help="bam file containing the reads")
parser.add_argument('--wd',type=str,default="contigAnalysis",help="the files will be  printed to this directory")
parser.add_argument('--cores',type=int, default =16, help="number of cores (default = 16)")
args= parser.parse_args()
main(args,wd)
