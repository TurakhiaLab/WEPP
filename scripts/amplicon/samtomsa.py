#!/usr/bin/python3 

# Generating contig file for reference
f=open("wuhCor1.fa","r")
a=f.readlines()
f.close()
ref=""
count=0
for i in a:
    read=i.split(" ")
    # print(read)
    if(count==0):
        refName=read[0].split(">")[1].split("\n")[0]
        count+=1
        continue
    else:
        ref+=read[0].split("\n")[0]
refLen=len(ref)

def decode(code):
    i=0
    length=list()
    type=list()
    while (True):
        if(code==""):
            break
        if(code[i]=="M" or code[i]=="S" or code[i]=="D" or code[i]=="I"):
            value=code[:i]
            length.append(value)
            type.append(code[i])
            code=code[i+1:]
            i=-1
        i+=1
    return (length, type)

fileread=open("data/amplicons.sam","r")
readLines=fileread.readlines()
fileread.close()

msafile=open("data/amplicons.fa","w")
msafile.write(">"+refName+"\n"+ref+"\n")
iterate=len(readLines)
for i in range(iterate):
    if(readLines[i][0]=="@"):
        continue
    else:
        
        column = readLines[i].split()
        #msafile.write(">"+column[0]+"\n")
        
        query=column[9]
        length, type=decode(column[5])
        refIndex=int(column[3])
        for j in range(len(type)):
            if(type[j]=="M"):
                break
            else:
                refIndex-=int(length[j])
        #print(column[0],refIndex)
        alignedSeq=""
        queryIndex=0
        initialInst=refIndex
        
        #print(length, type)
        alignedSeq+="-"*(initialInst-1)
        
        for j in range(len(type)):
            if(type[j]=="M"):
                # print(refIndex, int(length[j]))
                alignedSeq+=query[queryIndex:int(length[j])+queryIndex]
                queryIndex+=int(length[j])
            elif(type[j]=="S"):
                # print(refIndex, int(length[j]))
                alignedSeq+=query[queryIndex:int(length[j])+queryIndex]
                queryIndex+=int(length[j])
            elif(type[j]=="D"):
                alignedSeq+="-"*int(length[j])
                #queryIndex+=int(length[j])
            elif(type[j]=="I"):
#                #print(length[j])
#                alignedSeq+=queryIndex[queryIndex:int(length[j])+queryIndex]
                queryIndex+=int(length[j])
        # print(column[0],queryIndex)
        #msafile.write(">"+column[0]+" "+str(queryIndex)+" "+str(refIndex)+" "+ str(refIndex+queryIndex-1)+"\n")
        msafile.write(">"+column[0]+'\t'+str(queryIndex)+'\t'+str(refIndex)+'\t'+ str(refIndex+queryIndex-1)+"\n")
        alignedSeq+="-"*(refLen-len(alignedSeq))
        msafile.write(alignedSeq+"\n")
