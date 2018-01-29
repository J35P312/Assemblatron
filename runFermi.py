import os
import sys
import argparse

wd= os.path.dirname(os.path.realpath(__file__))
fermi="{}/fermikit/fermi.kit/".format(wd)

def main(args,fermi,wd):
    #create the fastq file if bam input 

    #apply bloom filter and build the index
    ropebwt="{}/ropebwt2 -m {} -dNCr - > {}.fmd 2> {}.fmd.log".format(fermi,args.batch,args.prefix,args.prefix)
    bfc="{}/bfc -1s {} -k {} -t {} {} 2> {}.flt.fq.gz.log".format(fermi,args.z,args.l,args.cores,args.fastq,args.prefix)
    os.system("{} | {}".format(bfc,ropebwt))
    #assemble and align
    os.system( "{}/fermi2 assemble -l {} -m {} -t {} {}.fmd 2> {}.pre.gz.log | gzip -1 > {}.pre.gz".format(fermi,args.l,args.m,args.cores,args.prefix,args.prefix,args.prefix) )
    os.system("{}/fermi2 simplify -CSo 66 -m {} -T 61 {}.pre.gz 2>  {}.mag.gz.log | bwa mem -X intractg -t {} {} - | samtools view -Sbh - | samtools sort - > {}.bam".format(fermi,args.m,args.prefix,args.prefix,args.cores,args.ref,args.prefix))
    print "{}/fermi2 simplify -CSo {} -m {} -T 61 {}.pre.gz 2>  {}.mag.gz.log | bwa mem -X intractg -t {} {} - | samtools view -Sbh - | samtools sort - > {}.bam".format(fermi,args.o,args.m,args.prefix,args.prefix,args.cores,args.ref,args.prefix) 
    os.system( "samtools index {}.bam".format(args.prefix) )

parser = argparse.ArgumentParser("""runFermi - a wrapper for the fermi assembler""")

parser.add_argument('--fastq',type=str, help="input fastq")
parser.add_argument('--ref',required = True,type=str, help="reference fasta")
parser.add_argument('--prefix',required = True,type=str, help="prefix of the output files")
parser.add_argument('--cores',type=int, default =16, help="number of cores (default = 16)")
parser.add_argument('--batch',type=str, default ="20g", help="batch size for multi-string indexing; 0 for single-string (default=20g)")
parser.add_argument('-z',type=str, default ="3G", help="genome size (use K,M, or G) (default = 3G)")
parser.add_argument('-l',type=int, default =51, help="min match (default = 51)")
parser.add_argument('-m',type=int, default =74, help="min merrge length (default = 74)")
parser.add_argument('-o',type=int, default =66, help="min merrge length (default = 74)")

args= parser.parse_args()
if args.fastq or args.bam:
    main(args,fermi,wd)
else:
    print "missing fastq or bam"
