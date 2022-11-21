import os
from Bio import SeqIO
from Bio.Seq import Seq
from ali import translate_nucseq
import sys

filediff=0
fileid=0

genelist_file =open("/home/mbastian/data/Zoonomia_postbusco/Genes/Genes_V2/genelist_filtred_V2","r")
genelist=genelist_file.readlines()
genelist = [line.rstrip("\n") for line in genelist]

spfile=open("/home/mbastian/data/Zoonomia_postbusco/Genes_to_analyse/180splist","r") #suposed the same as for the V2
splist = [line.rstrip("\n") for line in spfile]

output=open("/home/mbastian/data/Zoonomia_postbusco/seqidV1_V2toalign","w")
for gene in genelist:
    spdiff=0
    seqv1_dict = SeqIO.to_dict(SeqIO.parse("/home/mbastian/data/Zoonomia_postbusco/Genes_to_analyse/"+gene+"_filtred.fna", "fasta"))
    seqv2_dict = SeqIO.to_dict(SeqIO.parse("/home/mbastian/data/Zoonomia_postbusco/Genes_to_analyse_V2/" + gene + "_filtred.fna", "fasta"))
    listspv2=[]
    for sp2 in seqv2_dict :
        listspv2.append(sp2)
    for sp in seqv1_dict:
        if sp in listspv2:
            if str(seqv1_dict[sp].seq)!=str(seqv2_dict[sp].seq):
                #print("different seq")
                spdiff+=1
        else:
            spdiff+=1
    if spdiff !=0:
        filediff+=1



    else:
        fileid+=1
        output.write(gene + "\n")
print("nbr file diff = " + filediff)