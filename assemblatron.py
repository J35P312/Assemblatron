import os
import sys
import numpy
import argparse
#import readVCF
import itertools

def read_cigar(cigar):
    deletions=0
    insertions=0
    SC = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
    length=0
    first=True
    clip_after=True

    aligned_range=[]
    current_pos=1
    for i in range(0,len(SC)/2):
        if first and SC[i*2+1] == "M":
            first = False
        elif first and SC[i*2+1] == "S":
            first = False
            clip_after=False
        if SC[i*2+1] == "M":
            length += int( SC[i*2] )
            bases=range(0,int( SC[i*2] ))
            for j in range(0,len(bases)):
                bases[j] += current_pos

            aligned_range += bases
            current_pos += int( SC[i*2] )
       	elif SC[i*2+1] == "I":
       	    insertions+=1
            length += int( SC[i*2] )
            bases=range(0,int( SC[i*2] ))
            for j in range(0,len(bases)):
                bases[j] += current_pos
            aligned_range += bases

            current_pos += int( SC[i*2] )
       	elif SC[i*2+1] == "D":
       	    deletions +=1
        else:
            current_pos += int( SC[i*2] )

    return deletions,insertions,length,clip_after,aligned_range




def report_splits(bam,working_dir,args):
    headerString="##fileformat=VCFv4.1\n";
    headerString+="##source=ASSEMBLATRON\n";
    headerString+="##ALT=<ID=INV,Description=\"Inversion\">\n";
    headerString+="##ALT=<ID=BND,Description=\"Break end\">\n";
    headerString+="##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    headerString+="##INFO=<ID=END,Number=1,Type=String,Description=\"End of an intra-chromosomal variant\">\n";
    headerString+="##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
    print headerString
    
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


    bam_prefix=bam.split("/")[-1]
    prefix=working_dir + "/" +  bam_prefix[0:-4]
    calls={}
    print '{}_splits.sam'.format(prefix)
    for line in open( '{}_splits.sam'.format(prefix) ):

        content=line.strip().split("\t")
        contig_id=content[0]
        cigar=content[5]
        flag="{0:012b}".format(int(content[1]))
        if int(flag[-9]) or int(flag[0]):
            continue
        orientationA="+"
        if int(flag[-5]):
            orientationA="-"

        if args.len_ctg > len(content[9]):
            continue 
        cigar_del,cigar_ins,length,clip_after, range_primary=read_cigar(content[5])
        n_segments=0
        segments={}
        segments[n_segments] = {"contig_start_pos":min(range_primary),"contig_end_pos":max(range_primary),"contig_range":range_primary,"chr":content[2],"start": int(content[3]) ,"end": int(content[3])+length-1 ,"orientation":orientationA,"main":True}
        if orientationA == "-":
            segments[n_segments] = {"contig_start_pos":min(range_primary),"contig_end_pos":max(range_primary),"contig_range":range_primary,"chr":content[2],"end": int(content[3]) ,"start": int(content[3])+length-1 ,"orientation":orientationA,"main":True}
        if args.len > len(range_primary) or args.q > int(content[4]) :
            continue

        SA_line=line.strip().split("SA:Z:")[-1].split("\t")[0]
        SA_fields=SA_line.strip(";").split(";")
        for SA in SA_fields:
            
            split_read=SA.split(",")
            #print split_read
            pos = int(split_read[1])
            chromosome= split_read[0]

            cigar_del,cigar_ins,length,clip_after,range_secondary=read_cigar(split_read[3])
            if args.len > len(range_secondary) or args.q > int(split_read[4]):
                continue
            n_segments += 1
            SA_orientation=split_read[2]
            SA_start=pos
            SA_end=length+pos-1

            segments[n_segments] = {"contig_start_pos":min(range_secondary),"contig_end_pos":max(range_secondary),"contig_range":range_secondary,"chr":chromosome,"end": SA_end ,"start": SA_start ,"orientation":SA_orientation,"main":False}
            if SA_orientation == "-":
                segments[n_segments] = {"contig_start_pos":min(range_secondary),"contig_end_pos":max(range_secondary),"contig_range":range_secondary,"chr":chromosome,"end": SA_start ,"start": SA_end ,"orientation":SA_orientation,"main":False}


        if len(segments) == 1:
            continue
        distance_matrix=numpy.zeros( (len(segments),len(segments)))
        for segment in range(0,len(segments)):
            for neihbhour in range(0,len(segments)):
                distance_matrix[segment][neihbhour] = segments[neihbhour]["contig_start_pos"] - segments[segment]["contig_end_pos"]
                if(segment == neihbhour) or 0 > distance_matrix[segment][neihbhour] :
                    distance_matrix[segment][neihbhour] = float("inf")

        
        breakpoints={}
        for segment in range(0,len(segments)):
            closest_bp= min(distance_matrix[segment])
            if 0 > closest_bp or closest_bp == float("inf"):
                continue
            B=numpy.argmin(distance_matrix[segment])
            if not segments[segment]["chr"] in calls:
                calls[segments[segment]["chr"]]={}
            if not segments[B]["chr"] in calls[segments[segment]["chr"]]:
                calls[segments[segment]["chr"]][segments[B]["chr"]]={}
            calls[segments[segment]["chr"]][segments[B]["chr"]].add("{}\t{}\t{}\t{}\t{}\t{}".format(segments[segment]["chr"],segments[segment]["start"],segments[B]["chr"],segments[B]["start"],segments[segment]["orientation"],segments[B]["orientation"]))


    i=0
    for chrA in calls:
        for chrB in calls[chrA]:
            for var in calls[chrA][chrB]:
                content= var.split("\t")
                if chrA == chrB and content[-1] != content[-2]:
                    INFO="END={};SVTYPE={};SVLEN={}".format(content[3],"INV",abs(int(content[1])-int(content[3])))
                    vcf_line=[chrA,content[1],"var_{}".format(i),"N","<INV>",".","PASS",INFO]
                else:
                    if chrA == chrB:
                        INFO="SVTYPE={};SVLEN={}".format("BND",abs(int(content[1])-int(content[3])))
                    else:
                        INFO="SVTYPE={}".format("BND")
                    INFO+=";CONTIGID={};CIGAR={}".format(contig_id,cigar)
                    vcf_line=[chrA,content[1],"var_{}".format(i),"N","N[{}:{}[".format(content[2],content[3]),".","PASS",INFO]
                print "\t".join(vcf_line)
                i+=1

def find_splits(bam,working_dir,args):
    bam_prefix=bam.split("/")[-1]
    prefix=working_dir + "/" +  bam_prefix[0:-4]
    found = False
    #print "samtools view {} | grep \"SA:\" | grep -E :{},|;{}, > test.sam".format(bam)
    os.system("samtools view -h -q {} -F 256 {} | samtools view -hS -F 2048 - | samtools view -S -F 4 - | grep \"SA:\"  > {}_splits.sam".format(args.q,bam,prefix))
    if os.path.getsize( '{}_splits.sam'.format(prefix) ):
        found = True

    return found

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

