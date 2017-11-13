import os
import sys
import argparse

wd= os.path.dirname(os.path.realpath(__file__))
fermi="{}/fermikit/fermi.kit/".format(wd)

def main(args,fermi,wd):
    #build the index
    if args.bam:
          os.system("{} {} | {}/ropebwt2 -m {} -dNCr - > {}.fmd 2> {}.fmd.log".format(wd+"/bam2fq.sh",args.bam,fermi,args.batch,args.prefix,args.prefix))
    elif args.fastq:            
          os.system("{}/ropebwt2 -m {} -dNCr {} > {}.fmd 2> {}.fmd.log".format(fermi,args.fastq,args.batch,args.prefix,args.prefix))

    #assemble and align
    os.system( "{}/fermi2 assemble -l {} -m {} -t {} {}.fmd 2> {}.pre.gz.log | gzip -1 > {}.pre.gz".format(fermi,args.l,args.m,args.cores,args.prefix,args.prefix,args.prefix) )
    os.system("{}/fermi2 simplify -CSo 51 -m {} -T 51 {}.pre.gz 2>  {}.mag.gz.log | bwa mem -X intractg -t {} {} - | samtools view -Sbh - | samtools sort - {}".format(fermi,args.m,args.prefix,args.prefix,args.cores,args.ref,args.prefix))
    os.system( "samtools index {}.bam".format(args.prefix) )

parser = argparse.ArgumentParser("""runFermi - a wrapper for the fermi assembler""")

parser.add_argument('--fastq',type=str, help="input fastq")
parser.add_argument('--bam',type=str, help="input bam, may be used insead of fastq")
parser.add_argument('--ref',required = True,type=str, help="reference fasta")
parser.add_argument('--prefix',required = True,type=str, help="prefix of the output files")
parser.add_argument('--cores',type=int, default =16, help="number of cores (default = 16)")
parser.add_argument('--batch',type=str, default ="20g", help="batch size for multi-string indexing; 0 for single-string (default=20g)")
parser.add_argument('-l',type=int, default =41, help="min match (default = 41)")
parser.add_argument('-m',type=int, default =41, help="min merrge length (default = 41)")

args= parser.parse_args()
if args.fastq or args.bam:
    main(args,fermi,wd)
else:
    print "missing fastq or bam"
