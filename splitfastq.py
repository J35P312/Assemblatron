import sys
#accepts the number of output files as first argument, and an output prefix as second. The fastq input is accepted through piping
to_split=int(sys.argv[1])
file_handles=[]
prefix=sys.argv[2]
for i in range(0,to_split):
    file_handles.append(open("{}.{}.fastq".format(prefix,i),"w"))


itterator=0
file_to_print=0
for line in sys.stdin:
    file_handles[file_to_print].write(line)
    itterator+=1

    if itterator== 7:
        itterator=0
        file_to_print+=1
        if file_to_print == len(file_handles):
            file_to_print=0
