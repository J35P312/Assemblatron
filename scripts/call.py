import os
import numpy
import itertools
import math


def header():
    header=["##fileformat=VCFv4.1",
    "##source=Assemblatron",
    "##FILTER=<ID=MinLen,Description=\"Variant smaller than the minimum variant limit\">",
    "##FILTER=<ID=MaxCTG,Description=\"Too many contigs at the breakpoint\">",
    "##FILTER=<ID=MaxCov,Description=\"Too high coverage at the breakpoints\">",
    "##FILTER=<ID=ZeroPloidy,Description=\"less than 1X coverage across the chromosome\">",
    "##FILTER=<ID=MaxNeighbours,Description=\"too many calls withing a 1kb region\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    "##INFO=<ID=INSCTG,Number=2,Type=String,Description=\"Sequence of the insertion\">", 
    "##INFO=<ID=INSCIGAR,Number=1,Type=String,Description=\"Cigar operation of the insertion\">",
    "##INFO=<ID=CIGAR,Number=2,Type=String,Description=\"Cigar operation of the two alignments\">",
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
    "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Length of the inserted sequence\">",
    "##INFO=<ID=ORIENTATION,Number=2,Type=String,Description=\"Orientation of the two alignments\">",  
    "##INFO=<ID=MAPQ,Number=2,Type=Integer,Description=\"Mapping quality of the two alignments\">",
    "##INFO=<ID=ALNLEN,Number=2,Type=Integer,Description=\"Length(bp) of the two alignments\">", 
    "##INFO=<ID=CTGCOV,Number=2,Type=Float,Description=\"Number of contigs at the breakpoints\">",
    "##INFO=<ID=COV,Number=2,Type=Float,Description=\"coverage at the breakpoints\">",
    "##INFO=<ID=COVM,Number=1,Type=Float,Description=\"coverage between between the breakpoints\">",
    "##INFO=<ID=NEIGHBOURS,Number=2,Type=Float,Description=\"Number of SV contigs within a 500 bp radius\">",

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


def call_cigar_var(aln,min_len,max_ctg,coverage_structure):
    cigar_calls=[]
    length=0
    #compute the length of the contig
    pos=int(aln[3])
    aln_pos=0
    het_lim=1.5
    j=0
    cigar=aln[5]
    SC = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
    for i in range(0,len(SC)/2):
        if (SC[i*2+1] == "D" or SC[i*2+1] == "I")and int( SC[i*2] ) > min_len:
            length = int( SC[i*2] )
            cova=coverage_structure[aln[2]][int(math.floor(pos/100.0))]
            covb=coverage_structure[aln[2]][int(math.floor(pos/100.0))+int(math.floor(length/100.0))]

            if SC[i*2+1] == "I":
                ALT="<INS>"
                seq=aln[9][aln_pos:aln_pos+length]
                INFO="INSLEN={};INSCIGAR={};INSCTG={};SVTYPE=INS;MAPQ={},.;CIGAR={},.;CTGCOV={},{}".format(length,cigar,seq,aln[4],cigar,cova,covb);
                end=pos

            elif SC[i*2+1] == "D":
                ALT="<DEL>"
                INFO="END={};SVLEN={};SVTYPE=DEL;MAPQ={},.;CIGAR={},.;CTGCOV={},{}".format(length+pos-1,length,aln[4],cigar,cova,covb);
                end=length+pos-1

            FILTER="PASS"
            if covb >= max_ctg or cova >= max_ctg:
                FILTER="maxCTG"
            zygosity="1/1"
            if covb >= het_lim or cova >= het_lim:
                zygosity="0/1"

            VARID="{}_cigar_{}".format(aln[0],j)
            cigar_calls.append([pos , end, ALT ,"{}\t{}\t{}\tN\t{}\t.\t{}\t{}\tGT\t{}".format(aln[2],pos,VARID,ALT,FILTER,INFO,zygosity)])
            j+=1
        if (SC[i*2+1] == "D" or SC[i*2+1] == "M"):
            pos+=int( SC[i*2] )
        if (SC[i*2+1] == "I" or SC[i*2+1] == "M"):
            aln_pos+=int( SC[i*2] )

    return cigar_calls

def retrieve_var(coverage_structure,contigs,contig,var,order,i):

    chra=contigs[contig]["chr"][order[i]]

    if contigs[contig]["orientation"][order[i]] == "+":
        posa=contigs[contig]["pos"][order[i]]+compute_aln_length(contigs[contig]["cigar"][order[i]])-1
    else:
        posa=contigs[contig]["pos"][order[i]]

    if contigs[contig]["orientation"][order[i]] == "+":
        posb=contigs[contig]["pos"][order[i+1]]
    else:
        posb=contigs[contig]["pos"][order[i+1]]+compute_aln_length(contigs[contig]["cigar"][order[i+1]])-1
    
    index_a=int(math.floor(posa/100.0))
    index_b=int(math.floor(posb/100.0))
    if index_a-1 > 0:
        cova=max(coverage_structure[contigs[contig]["chr"][order[i]]] [ index_a-1:index_a+2] )
    else:
        cova=max(coverage_structure[contigs[contig]["chr"][order[i]]] [ 0:index_a+2] )
    if index_b-1 > 0:
        covb=max(coverage_structure[contigs[contig]["chr"][order[i+1]]] [ index_b-1:index_b+2] )
    else:
        covb=max(coverage_structure[contigs[contig]["chr"][order[i+1]]] [ 0:index_b+2] )

    return({"type":var,"chrA":chra,"chrB":contigs[contig]["chr"][order[i+1]],"start":posa,"end":posb,"oa":contigs[contig]["orientation"][order[i]],"ob":contigs[contig]["orientation"][order[i+1]],"contig":contig,"qa":contigs[contig]["q"][order[i]],"qb":contigs[contig]["q"][order[i+1]],"cigara":contigs[contig]["cigar"][order[i]],"cigarb":contigs[contig]["cigar"][order[i+1]],"call":i,"lena":compute_aln_length(contigs[contig]["cigar"][order[i]]),"lenb":compute_aln_length(contigs[contig]["cigar"][order[i+1]]),"cova":cova , "covb":covb })

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

    ascii= {"\"":1,"#":2,"$":3,"%":4,"&":5,"\'":6,"(":7,")":8,"*":9,"+" :10,"," :11,"-" :12,".":13,"/":14,"0":15,"1":16,"2":17,"3":18,"4":19,"5":20,"6":21,"7":22,"8":23,"9":24,":":25,";":26,"<":27,"=":28,">":29,"?":30,"@":31,"A":32,"B":33,"C":34,"D":35,"E":36,"F":37,"G":38,"H":39,"I":40,"J":41,"K":42,"L":43,"M":44,"N":45,"O":46,"P":47,"Q":48,"R":49,"S":50,"T":51,"U":52,"V":53,"W":54,"X":55,"Y":56,"Z":57,"[":58,"\\":59,"]":60,"^":61,"_":62,"`":63,"a":64,"b":65,"c":66,"d":67,"e":68,"f":69,"g":70,"h":71,"i":72,"j":73,"k":74,"l":75,"m":76,"n":77,"o":78,"p":79,"q":80,"r":81,"s":82,"t":83,"u":84,"q":85,"r":86,"s":87,"t":88,"u":89,"v":90,"w":91,"x":92,"y":93,"z":94,"{":95,"|":96,"}":97,"~":98}

    ploidy=2
    kmer=args.len_ctg
    q=args.q
    min_len=args.min_size
    complex_dist=5000
    ins_dist=100
    duplicate_dist=10
    max_ctg=args.max_contigs

    if args.sample:
        sample=args.sample
    else:
        sample=args.bam.split("/")[-1].split(".")[0]

    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(sample)

    coverage_structure={}
    coverage_contig_coverage={}
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
                    coverage_contig_coverage[chromosome]=numpy.zeros(bins)
                    chromosome_order.append(chromosome)

    contigs={}
    variant_calls={}
    with os.popen("samtools view {}".format(args.bam)) as sam:
        for line in sam:
            content=line.strip().split()

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
                        filled_elements=int(math.floor(sequence_length/100.0))
                        remainingRead=(sequence_length % binSize)/100
                        contig_pos=0
                        for i in range(0,filled_elements):
                            coverage_structure[content[2]][element]+=1;
                            
                            if contig_pos+100 >= len(content[10]):
                                coverage_contig_coverage[content[2]][element]+= sum( ascii[idx] for idx in content[10][contig_pos:len(content[10])] ) /100.0
                            else:
                                coverage_contig_coverage[content[2]][element]+= sum( ascii[idx] for idx in content[10][contig_pos:contig_pos+100] ) /100.0
                            contig_pos+=100
                            element+=1

                        if (remainingRead > 0 and len(coverage_structure[content[2]]) > element+1):
                            coverage_structure[content[2]][element]+=remainingRead/100.0;
                            coverage_contig_coverage[content[2]][element]+= sum( ascii[idx] for idx in content[10][contig_pos:len(content[10])] ) /100.0

            if not "\tSA:Z" in line:
                continue

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
    
    mean_coverage=0
    n_bins=0
    chrom_cov={}
    ploidies={}
    for chromosome in coverage_contig_coverage:
            idx=numpy.where( coverage_contig_coverage[chromosome] > 0 )
            chrom_bins=sum(coverage_contig_coverage[chromosome][idx])
            n_bins_chrom=len(idx[0])
            chrom_cov[chromosome]=numpy.median(coverage_contig_coverage[chromosome])
            mean_coverage+=chrom_bins
            n_bins+= n_bins_chrom

    mean_coverage=mean_coverage/float(n_bins)
    for chromosome in coverage_contig_coverage:
        ploidies[chromosome]=ploidy*chrom_cov[chromosome]/mean_coverage
        #print [ploidies[chromosome],chrom_cov[chromosome],chromosome,mean_coverage]

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
                        if "+" == contigs[contig]["orientation"][i]:
                            connectivity[k][0]+=1
                            connectivity[k][i+1]+=1
                        else:
                            connectivity[contig_length-k-1][0]+=1
                            connectivity[contig_length-k-1][i+1]+=1
                    
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
               
                if len(order) == 3 and contigs[contig]["q"][order[0]] >= q and contigs[contig]["q"][order[-1]] >= q and contigs[contig]["chr"][order[0]] == contigs[contig]["chr"][order[-1]] and abs(contigs[contig]["pos"][order[0]] - contigs[contig]["pos"][order[-1]]) < ins_dist and (contigs[contig]["chr"][order[0]] != contigs[contig]["chr"][order[1]] or abs(contigs[contig]["pos"][order[0]]-contigs[contig]["pos"][order[1]]) > ins_dist ):
                    calls.append(retrieve_var(coverage_structure,contigs,contig,"INS",order,i))
                    calls[-1]["oc"]=contigs[contig]["orientation"][order[-1]]
                    calls[-1]["qc"]=contigs[contig]["q"][order[-1]]
                    calls[-1]["cigarc"]=contigs[contig]["cigar"][order[-1]]
                    calls[-1]["lenc"]=compute_aln_length(contigs[contig]["cigar"][order[-1]]) 
                    calls[-1]["seqb"]=contigs[contig]["sequence"][order[1]]
                    idx_a=int(math.floor(calls[-1]["start"]/100.0))
                    idx_b=int(math.floor(calls[-1]["end"]/100.0))
                    calls[-1]["covA"]=coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_a ]
                    calls[-1]["covB"]=coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_a ]
                    break

                 #discard low quality alignments
                if contigs[contig]["q"][order[i]] < q or contigs[contig]["q"][order[i+1]] < q:
                    continue

                if not contigs[contig]["chr"][order[i]] == contigs[contig]["chr"][order[i+1]]:
                    calls.append(retrieve_var(coverage_structure,contigs,contig,"BND",order,i))
                    idx_a=int(math.floor(calls[-1]["start"]/100.0))
                    idx_b=int(math.floor(calls[-1]["end"]/100.0))
                    calls[-1]["covA"]=coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_a ]
                    calls[-1]["covB"]=coverage_contig_coverage[ contigs[contig]["chr"][order[i+1]] ][ idx_b ]

                elif contigs[contig]["orientation"][order[i]] == contigs[contig]["orientation"][order[i+1]]:

                    calls.append(retrieve_var(coverage_structure,contigs,contig,"BND",order,i)) 

                    if calls[-1]["end"] < calls[-1]["start"]:
                        tmp=calls[-1]["end"]
                        calls[-1]["end"]=calls[-1]["start"]
                        calls[-1]["start"]=tmp


                    idx_a=int(math.floor(calls[-1]["start"]/100.0))
                    idx_b=int(math.floor(calls[-1]["end"]/100.0))

                    cov_between=numpy.average(coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_a:idx_b+1 ])
                    calls[-1]["cov_between"]=cov_between
                    calls[-1]["covA"]=coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_a ]
                    calls[-1]["covB"]=coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_b ]
                   
                    gain=False
                    loss=False
                    if ploidies[ contigs[contig]["chr"][order[i]] ] > 0:
                        if cov_between/chrom_cov[ contigs[contig]["chr"][order[i]] ] <  abs(1-0.5/float(ploidies[ contigs[contig]["chr"][order[i]] ])):
                            loss=True
                        elif cov_between/chrom_cov[ contigs[contig]["chr"][order[i]] ] >  1+0.5/float(ploidies[ contigs[contig]["chr"][order[i]] ]):
                            gain=True
                    
                        if contigs[contig]["orientation"][order[i]] == "+":
                            if contigs[contig]["pos"][order[i]] < contigs[contig]["pos"][order[i+1]]:
                                if loss:
                                    calls[-1]["type"]="DEL"
                            else:
                                if gain:
                                    calls[-1]["type"]="TDUP"
                        else:
                            if contigs[contig]["pos"][order[i]] > contigs[contig]["pos"][order[i+1]]:
                                if loss:
                                    calls[-1]["type"]="DEL"
                            else:
                                if gain:
                                    calls[-1]["type"]="TDUP"
                else:

                    calls.append(retrieve_var(coverage_structure,contigs,contig,"INV",order,i))

                    if calls[-1]["end"] < calls[-1]["start"]:
                        tmp=calls[-1]["end"]
                        calls[-1]["end"]=calls[-1]["start"]
                        calls[-1]["start"]=tmp

                    idx_a=int(math.floor(calls[-1]["start"]/100.0))
                    idx_b=int(math.floor(calls[-1]["end"]/100.0))
                    cov_between=numpy.average(coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_a:idx_b+1 ])
                    calls[-1]["cov_between"]=cov_between
                    calls[-1]["covA"]=coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_a ]
                    calls[-1]["covB"]=coverage_contig_coverage[ contigs[contig]["chr"][order[i]] ][ idx_b ]



    breakpoints={}
    contig_ids={}
    for call in calls:
        if not call["chrA"] in breakpoints:
            breakpoints[call["chrA"]]=[]
            contig_ids[call["chrA"]]=[]            

        if not call["chrB"] in breakpoints:
            breakpoints[call["chrB"]]=[]
            contig_ids[call["chrB"]]=[]

                    
        breakpoints[call["chrA"]].append(call["start"])
        breakpoints[call["chrB"]].append(call["end"])

        contig_ids[call["chrA"]].append(call["contig"])
        contig_ids[call["chrB"]].append(call["contig"])


        
    for chromosome in breakpoints:
        breakpoints[chromosome]=numpy.array(breakpoints[chromosome])
        contig_ids[chromosome]=numpy.array(contig_ids[chromosome])

    het_lim=1.5
    for call in calls:
        if not call["chrA"] in variant_calls:
            variant_calls[call["chrA"]]=[]

        if call["type"] == "BND":
            bps_A=numpy.where( (breakpoints[call["chrA"]] > call["start"]-500) & (breakpoints[call["chrA"]] < call["start"]+500)  )
            bps_B=numpy.where( (breakpoints[call["chrB"]] > call["end"]-500) & (breakpoints[call["chrB"]] < call["end"]+500)  )
    
            contig_ids_a=set(contig_ids[call["chrA"]][bps_A])
            contig_ids_b=set(contig_ids[call["chrB"]][bps_B])

            if args.skip_inter and call["chrA"] != call["chrB"]:
                continue


            INFO="SVTYPE={};MAPQ={},{};CIGAR={},{};ORIENTATION={},{},ALNLEN={},{};CTGCOV={},{},;COV={},{};NEIGHBOURS={},{}".format(call["type"],call["qa"],call["qb"],call["cigara"],call["cigarb"],call["oa"],call["ob"],call["lena"],call["lenb"],call["cova"],call["covb"],call["covA"],call["covB"],len(contig_ids_a),len(contig_ids_b));
            VARID="{}_{}".format(call["contig"],call["call"])
        
            if call["ob"] == call["ob"]:
                ALT="N[{}:{}[".format(call["chrB"],call["end"])
            else:
                ALT="[{}:{}[N".format(call["chrB"],call["end"])
    
            FILTER="PASS"
            if call["covb"] >= max_ctg or call["cova"] >= max_ctg:
                FILTER="maxCTG"
            if call["covA"] >= args.max_coverage*chrom_cov[ call["chrA"] ] or call["covB"] >= args.max_coverage*chrom_cov[ call["chrB"] ]:
                FILTER="MaxCov"
            if  1 > chrom_cov[ call["chrA"] ] or  1 > chrom_cov[ call["chrB"] ]:
                FILTER="ZeroPloidy"
            if len(contig_ids_a) > 3 or len(contig_ids_b) > 3:
                FILTER="MaxNeighbours"

            zygosity="1/1"
            if call["covb"] >= het_lim or call["cova"] >= het_lim:
                zygosity="0/1"

            variant_calls[call["chrA"]].append([call["start"],call["end"],"BND","{}\t{}\t{}\tN\t{}\t.\t{}\t{}\tGT\t{}".format(call["chrA"],call["start"],VARID,ALT,FILTER,INFO,zygosity)])
        elif call["type"] == "INS":
            bps_A=numpy.where( (breakpoints[call["chrA"]] > call["start"]-500) & (breakpoints[call["chrA"]] < call["start"]+500)  )
            bps_B=numpy.where( (breakpoints[call["chrB"]] > call["end"]-500) & (breakpoints[call["chrB"]] < call["end"]+500)  )
    
            contig_ids_a=set(contig_ids[call["chrA"]][bps_A])
            contig_ids_b=set(contig_ids[call["chrB"]][bps_B])

            INFO="INSLEN={};INSCIGAR={};INSCTG={};SVTYPE={};MAPQ={},{};CIGAR={},{};ORIENTATION={},{},ALNLEN={},{};CTGCOV={},{};COV={},{};NEIGHBOURS={},{}".format(call["lenb"],call["cigarb"],call["seqb"],call["type"],call["qa"],call["qc"],call["cigara"],call["cigarc"],call["oa"],call["oc"],call["lena"],call["lenc"],call["cova"],call["cova"],call["covA"],call["covA"],len(contig_ids_a),len(contig_ids_b));

            VARID="{}_{}".format(call["contig"],call["call"])
            ALT="<{}>".format(call["type"])
            FILTER="PASS"
            length= abs(call["start"]-call["end"])+1
            if length < min_len:
                FILTER="MinLen"
            if call["cova"] >= max_ctg:
                FILTER="maxCTG"
            if call["covA"] >= args.max_coverage*chrom_cov[ call["chrA"] ]:
                FILTER="MaxCov"
            if  1 > chrom_cov[ call["chrA"]]:
                FILTER="ZeroPloidy"
            if len(contig_ids_a) > 3 or len(contig_ids_b) > 3:
                FILTER="MaxNeighbours"

            zygosity="1/1"
            if call["covb"] >= het_lim or call["cova"] >= het_lim:
                zygosity="0/1"
            variant_calls[call["chrA"]].append([call["start"],call["end"],ALT,"{}\t{}\t{}\tN\t{}\t.\t{}\t{}\tGT\t{}".format(call["chrA"],call["start"],VARID,ALT,FILTER,INFO,zygosity)])

        else: 
            bps_A=numpy.where( (breakpoints[call["chrA"]] > call["start"]-500) & (breakpoints[call["chrA"]] < call["start"]+500)  )
            bps_B=numpy.where( (breakpoints[call["chrB"]] > call["end"]-500) & (breakpoints[call["chrB"]] < call["end"]+500)  )
            contig_ids_a=set(contig_ids[call["chrA"]][bps_A])
            contig_ids_b=set(contig_ids[call["chrB"]][bps_B])

            INFO="END={};SVLEN={};SVTYPE={};MAPQ={},{};CIGAR={},{};ORIENTATION={},{},ALNLEN={},{};CTGCOV={},{};COV={},{};COVM={};NEIGHBOURS={},{}".format(call["end"],abs(call["start"]-call["end"])+1,call["type"],call["qa"],call["qb"],call["cigara"],call["cigarb"],call["oa"],call["ob"],call["lena"],call["lenb"],call["cova"],call["covb"],call["covA"],call["covB"],call["cov_between"],len(contig_ids_a),len(contig_ids_b));
            VARID="{}_{}".format(call["contig"],call["call"])
            ALT="<{}>".format(call["type"])
            FILTER="PASS"
            length= abs(call["start"]-call["end"])+1
            if length < min_len:
                FILTER="MinLen"
            if call["covb"] >= max_ctg or call["cova"] >= max_ctg:
                FILTER="maxCTG"

            if call["covA"] >= args.max_coverage*chrom_cov[ call["chrA"] ] or call["covB"] >= args.max_coverage*chrom_cov[ call["chrB"] ]:
                FILTER="MaxCov"
            if  1 > chrom_cov[ call["chrA"]] or  1 > chrom_cov[ call["chrB"] ]:
                FILTER="ZeroPloidy"
            if len(contig_ids_a) > 3 or len(contig_ids_b) > 3:
                FILTER="MaxNeighbours"

            zygosity="1/1"
            if call["covb"] >= het_lim or call["cova"] >= het_lim:
                zygosity="0/1"

            variant_calls[call["chrA"]].append([call["start"],call["end"],ALT,"{}\t{}\t{}\tN\t{}\t.\t{}\t{}\tGT\t{}".format(call["chrA"],call["start"],VARID,ALT,FILTER,INFO,zygosity)])


    with os.popen("samtools view {} -m 3 -q {}".format(args.bam,q)) as sam:
        for line in sam:
            content=line.strip().split()

            if not content[2] in variant_calls:
                variant_calls[content[2]]=[]

            call_list=call_cigar_var(content,min_len,max_ctg,coverage_structure)   
            variant_calls[content[2]]+=call_list 
                
    for chromosome in chromosome_order:
        first=True
        if not chromosome in variant_calls:
            continue
        for call in sorted(variant_calls[chromosome],key=lambda x: x[0]):
            if first:
                last_call=call
                print call[-1]
                first=False
                continue

            if abs(call[0]-last_call[0]) < duplicate_dist and abs(call[1]-last_call[1]) < duplicate_dist and call[2] == last_call[2]:
                continue

            else: 
                last_call=call
                print call[-1]

