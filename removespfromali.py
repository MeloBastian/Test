import sys
import ali

genelist_file =open(sys.argv[1],"r")
genelist=genelist_file.readlines()
genelist = [line.rstrip("\n") for line in genelist]

spfile=open("/beegfs/data/mbastian/Zoonomia_postbusco/Genes_to_analyse_V2/180sp","r")
splist = [line.rstrip("\n") for line in spfile]

for gene in genelist:
    aliset = ali.SequenceAlignment(file_name="/beegfs/data/mbastian/Zoonomia_postbusco/Genes/Genes_V2/"+gene+"_no_filter/"+ gene +".fasta", format="fasta")
    aliset.redefine_taxon_set(splist, all_taxa=False)
    aliset.write_fasta_to_file("/beegfs/data/mbastian/Zoonomia_postbusco/Genes_to_analyse_V2/"+gene+"_filtred.fna")




