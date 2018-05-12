import os
import numpy
import itertools
import math

def header():
    header=["##fileformat=VCFv4.1",
    "##source=Assemblatron",
    "##FILTER=<ID=MinLen,Description=\"Variant smaller than the minimum variant limit\">",
    "##FILTER=<ID=MaxCov,Description=\"Coverage higher than the coverage limit\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    "##INFO=<ID=INSCTG,Number=2,Type=String,Description=\"Sequence of the insertion\">", 
    "##INFO=<ID=INSCIGAR,Number=1,Type=String,Description=\"Cigar operation of the insertion\">",
    "##INFO=<ID=CIGAR,Number=2,Type=String,Description=\"Cigar operation of the two alignments\">",
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the inserted sequence\">",
    "##INFO=<ID=ORIENTATION,Number=2,Type=String,Description=\"Orientation of the two alignments\">",  
    "##INFO=<ID=MAPQ,Number=2,Type=Integer,Description=\"Mapping quality of the two alignments\">",
    "##INFO=<ID=ALNLEN,Number=2,Type=Integer,Description=\"Length(bp) of the two alignments\">", 
    "##INFO=<ID=COV,Number=2,Type=Float,Description=\"coverage at the breakpoints\">", 

]
    print "\n".join(header)

def compute_aln_length(cigar):
    length=0
    #compute the length of the contig
    SC = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
    for i in range(0,len(SC)/2):
        if SC[i*2+1] == "M":
            length += int( SC[i*2] )
    return length

def retrieve_var(coverage_structure,contigs,contig,var,order,i):

    chra=contigs[contig]["chr"][order[i]]
    posa=contigs[contig]["pos"][order[i]]+compute_aln_length(contigs[contig]["cigar"][order[i]])-1
    index_a=int(math.floor(contigs[contig]["pos"][order[i]]/100.0))+int(math.floor(compute_aln_length(contigs[contig]["cigar"][order[i]])/100.0))

    return({"type":var,"chrA":chra,"chrB":contigs[contig]["chr"][order[i+1]],"start":posa,"end":contigs[contig]["pos"][order[i+1]],"oa":contigs[contig]["orientation"][order[i]],"ob":contigs[contig]["orientation"][order[i+1]],"contig":contig,"qa":contigs[contig]["q"][order[i]],"qb":contigs[contig]["q"][order[i+1]],"cigara":contigs[contig]["cigar"][order[i]],"cigarb":contigs[contig]["cigar"][order[i+1]],"call":i,"lena":compute_aln_length(contigs[contig]["cigar"][order[i]]),"lenb":compute_aln_length(contigs[contig]["cigar"][order[i+1]]),"cova":coverage_structure[chra][ index_a ] , "covb":coverage_structure[contigs[contig]["chr"][order[i+1]]] [ int(math.floor(contigs[contig]["pos"][order[i+1]]/100.0)) ] })

def order_segments(assigned_segments,kmer):
        #remove segments that are too short, and determine the order of segments
        order=[]
        first=True
        for i in range(0,len(assigned_segments)):
            if first:
                current_segment=assigned_segments[i]
                length=1
                segment_start_pos=0
                first=False
                continue

            if current_segment == assigned_segments[i] and i < len(assigned_segments) -1:
                length+=1
            elif current_segment == assigned_segments[i] and current_segment !=-1 and i == len(assigned_segments) -1:
                    if not order and length > kmer:
                        order.append(current_segment)
                    elif not order:
                        pass
                    elif order[-1] != current_segment and length > kmer:
                        order.append(current_segment)
            else:
                if current_segment != -1:
                    if not order and length > kmer:
                        order.append(current_segment)
                    elif not order:
                        pass
                    elif order[-1] != current_segment and length > kmer:
                        order.append(current_segment)

                current_segment=assigned_segments[i]
                segment_start_pos=i
                length=1

	return(order)


def main(args):
    header()

    kmer=args.len_ctg
    q=args.q
    min_len=args.min_size
    complex_dist=5000
    ins_dist=100
    max_cov=args.max_coverage
    if args.sample:
        sample=args.sample
    else:
        sample=args.bam.split("/")[-1].split(".")[0]
    contigs={}


    coverage_structure={}
    chromosome_order=[]
    with os.popen("samtools view -H {}".format(args.bam)) as pipe:
        for line in pipe:
            if line[0] == "@":
                if "SN:" in line:
                    content=line.strip().split()
                    chromosome=content[1].split("SN:")[-1]
                    length=int(content[2].split("LN:")[-1])
                    bins=int( math.ceil(length/float(100)) )
                    coverage_structure[chromosome]=numpy.zeros(bins)
                    chromosome_order.append(chromosome)

    with os.popen("samtools view {} | grep SA:Z".format(args.bam)) as sam:
        for line in sam:
            content=line.strip().split()
    
            if not content[0] in contigs:
                contigs[content[0]]={}
                contigs[content[0]]["pos"]=[]
                contigs[content[0]]["chr"]=[]
                contigs[content[0]]["orientation"]=[]
                contigs[content[0]]["cigar"]=[]
                contigs[content[0]]["q"]=[]
                contigs[content[0]]["sequence"]=[]

            contigs[content[0]]["pos"].append(int(content[3]))
            contigs[content[0]]["chr"].append(content[2])
    
            orientation="+"
    
            flag="{0:012b}".format(int(content[1]))
            if int(flag[-5]):
                orientation= "-"
    
            contigs[content[0]]["orientation"].append(orientation)
            contigs[content[0]]["cigar"].append(content[5])
            contigs[content[0]]["q"].append(int(content[4]))
            contigs[content[0]]["sequence"].append(content[9])
    
            if content[2] in chromosome_order:
                sequence_length=compute_aln_length(content[5])
                binSize=100
                element=int(math.floor(int(content[3])/100.0))
                

                #coverage_structure[content[2]][pos]

                #if the entire read is inside the region, add all the bases to sequenced bases
                if(int(content[3]) >= element*binSize and int(content[3])+sequence_length-1 <= (element+1)*binSize):
                        coverage_structure[content[2]][element] += sequence_length/100.0;
                else:
                #if the read starts within the region but reaches outside it, add only those bases that fit inside the region.
                        coverage_structure[content[2]][element]+=((element+1)*binSize-int(content[3])+1)/float(binSize)

                        #the part of the read hanging out of the bin is added to the bins following the currentbin
                        remainingRead=sequence_length-((element+1)*binSize-int(content[3])+1);
                        while (remainingRead >= binSize and  len(coverage_structure[content[2]]) > element+1 ):
                                element+=1;
                                coverage_structure[content[2]][element]+=1;
                                remainingRead=remainingRead-binSize;
                        
                        if (remainingRead > 0 and len(coverage_structure[content[2]]) > element+1):
                                element+=1;
                                coverage_structure[content[2]][element]+=remainingRead/100.0;
                        





    #create alignment "maps" across all the contigs
    calls=[]
    for contig in contigs:
        splits=len(contigs[contig]["pos"])
    
        if splits < 2:
            continue

        contig_length=0
        #compute the length of the contig
        SC = ["".join(x) for _, x in itertools.groupby(contigs[contig]["cigar"][0], key=str.isdigit)]
        for i in range(0,len(SC)/2):
            if SC[i*2+1] == "M" or SC[i*2+1] == "S" or SC[i*2+1] == "H" or SC[i*2+1] == "I":
                contig_length += int( SC[i*2] )

        #create a table that describes where the alignments are located on the contig
        connectivity=numpy.zeros((contig_length,splits+1))
        for i in range(0,splits):
            SC = ["".join(x) for _, x in itertools.groupby(contigs[contig]["cigar"][i], key=str.isdigit)]
            pos=0
            for j in range(0,len(SC)/2):
            
                if SC[j*2+1] == "M" or SC[j*2+1] == "I":
                    length = int( SC[j*2] )
                    #print [contig,length,pos,contigs[contig]["cigar"][i]]
                    for k in range(pos,pos+length):
                        connectivity[k][0]+=1
                        connectivity[k][i+1]+=1
                if SC[j*2+1] != "D":
                    pos+=int( SC[j*2] )

        segments=[]
        #find regions on the contig that are aligned uniquely
        for i in range(0,contig_length):
            if connectivity[i][0] == 1:
                segments.append(1)
            else:
                segments.append(0)

        # now assign the unique regions to the alignments
        assigned_segments=[]
        for i in range(0,contig_length):
            if segments[i]:
                for j in range(0,splits):
                    if connectivity[i][j+1]:
                        assigned_segments.append(j)
            else:
                assigned_segments.append(-1)

        order=order_segments(assigned_segments,kmer)
    
        if len(order) > 1:
    
            for i in range(0,len(order)-1):
                #discard low quality alignments
                if len(order) == 3 and contigs[contig]["q"][order[0]] >= q and contigs[contig]["q"][order[-1]] >= q and contigs[contig]["chr"][order[0]] == contigs[contig]["chr"][order[-1]] and abs(contigs[contig]["pos"][order[0]] - contigs[contig]["pos"][order[-1]]) < ins_dist and (contigs[contig]["chr"][order[0]] != contigs[contig]["chr"][order[1]] or abs(contigs[contig]["pos"][order[0]]-contigs[contig]["pos"][order[1]]) > ins_dist ):
                    calls.append(retrieve_var(coverage_structure,contigs,contig,"INS",order,i))
                    calls[-1]["oc"]=contigs[contig]["orientation"][order[-1]]
                    calls[-1]["qc"]=contigs[contig]["q"][order[-1]]
                    calls[-1]["cigarc"]=contigs[contig]["cigar"][order[-1]]
                    calls[-1]["lenc"]=compute_aln_length(contigs[contig]["cigar"][order[-1]]) 
                    calls[-1]["seqb"]=contigs[contig]["sequence"][order[1]]
                    break

                if contigs[contig]["q"][order[i]] < q or contigs[contig]["q"][order[i+1]] < q:
                    continue

                if contig_length/float(len(order)) < complex_dist and len(order) > 2:
                    calls.append(retrieve_var(coverage_structure,contigs,contig,"BND",order,i))                
                elif not contigs[contig]["chr"][order[i]] == contigs[contig]["chr"][order[i+1]]:
                    calls.append(retrieve_var(coverage_structure,contigs,contig,"BND",order,i))
                elif contigs[contig]["orientation"][i] == contigs[contig]["orientation"][order[i+1]]:
                    if contigs[contig]["pos"][order[i]] < contigs[contig]["pos"][order[i+1]]:
                        calls.append(retrieve_var(coverage_structure,contigs,contig,"DEL",order,i))
                    else:
                        calls.append(retrieve_var(coverage_structure,contigs,contig,"TDUP",order,i))
                else:
                    calls.append(retrieve_var(coverage_structure,contigs,contig,"INV",order,i))

    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(sample)
    het_lim=1.5
    for call in calls:
        if call["type"] == "BND":
            INFO="SVTYPE={};MAPQ={},{};CIGAR={},{};ORIENTATION={},{},ALNLEN={},{};COV={},{}".format(call["type"],call["qa"],call["qb"],call["cigara"],call["cigarb"],call["oa"],call["ob"],call["lena"],call["lenb"],call["cova"],call["covb"]);
            VARID="{}_{}".format(call["contig"],call["call"])
        
            if call["ob"] == call["ob"]:
                ALT="N[{}:{}[".format(call["chrB"],call["end"])
            else:
                ALT="[{}:{}[N".format(call["chrB"],call["end"])
    
            FILTER="PASS"
            if call["covb"] >= max_cov or call["cova"] >= max_cov:
                FILTER="MaxCov"
            zygosity="1/1"
            if call["covb"] >= het_lim or call["cova"] >= het_lim:
                zygosity="0/1"
            print "{}\t{}\t{}\tN\t{}\t.\t{}\t{}\tGT\t{}".format(call["chrA"],call["start"],VARID,ALT,FILTER,INFO,zygosity)
        elif call["type"] == "INS":
            INFO="INSLEN={};INSCIGAR={};INSCTG={};SVTYPE={};MAPQ={},{};CIGAR={},{};ORIENTATION={},{},ALNLEN={},{};COV={},{}".format(call["lenb"],call["cigarb"],call["seqb"],call["type"],call["qa"],call["qc"],call["cigara"],call["cigarc"],call["oa"],call["oc"],call["lena"],call["lenc"],call["cova"],call["cova"]);

            VARID="{}_{}".format(call["contig"],call["call"])
            ALT="<{}>".format(call["type"])
            FILTER="PASS"
            length= abs(call["start"]-call["end"])+1
            if length < min_len:
                FILTER="MinLen"
            if call["cova"] >= max_cov:
                FILTER="MaxCov"
            zygosity="1/1"
            if call["covb"] >= het_lim or call["cova"] >= het_lim:
                zygosity="0/1"
            print "{}\t{}\t{}\tN\t{}\t.\t{}\t{}\tGT\t{}".format(call["chrA"],call["start"],VARID,ALT,FILTER,INFO,zygosity)

        else: 
            INFO="END={};SVLEN={};SVTYPE={};MAPQ={},{};CIGAR={},{};ORIENTATION={},{},ALNLEN={},{};COV={},{}".format(call["end"],abs(call["start"]-call["end"])+1,call["type"],call["qa"],call["qb"],call["cigara"],call["cigarb"],call["oa"],call["ob"],call["lena"],call["lenb"],call["cova"],call["covb"]);
            VARID="{}_{}".format(call["contig"],call["call"])
            ALT="<{}>".format(call["type"])
            FILTER="PASS"
            length= abs(call["start"]-call["end"])+1
            if length < min_len:
                FILTER="MinLen"
            if call["covb"] >= max_cov or call["cova"] >= max_cov:
                FILTER="MaxCov"
            zygosity="1/1"
            if call["covb"] >= het_lim or call["cova"] >= het_lim:
                zygosity="0/1"
            print "{}\t{}\t{}\tN\t{}\t.\t{}\t{}\tGT\t{}".format(call["chrA"],call["start"],VARID,ALT,FILTER,INFO,zygosity)
    
