import sys
import os
import ali


aliset = ali.MultiGeneSequenceAlignment(list_name=sys.argv[1], dir_name=os.getcwd(), format="fasta")

print("before filtering")
for gene in aliset.get_gene_list():
    print(gene, aliset.get_gene_ali(gene).get_ntaxa(), aliset.get_gene_ali(gene).get_nsite())

#aliset.prune_min_ntaxa(42)
aliset.prune_min_nsite(500)

print()
print("after filtering")
for gene in aliset.get_gene_list():
    print(gene, aliset.get_gene_ali(gene).get_ntaxa(), aliset.get_gene_ali(gene).get_nsite())

print()
taxon2nbrgene=dict()
for taxon in aliset.taxon_set:
        taxon2nbrgene[taxon]=0
for gene, value in aliset.alignments.items():
        for key in taxon2nbrgene.keys():
                if key in value.get_taxon_set():
                        taxon2nbrgene[key]+=1
print(taxon2nbrgene)

print()
aliset.homogeneize_taxon_sets()

print()
print("after homogeneization")
for gene in aliset.get_gene_list():
    print(gene, aliset.get_gene_ali(gene).get_ntaxa(), aliset.get_gene_ali(gene).get_nsite())

print()

print(len(aliset.alignments))

print("write all alignments to file")
aliset.write_all_to_file(sys.argv[2])
print("concatenate")
concat = aliset.concatenate()
print("total number of sites: ", concat.get_nsite())
print("write concat to file")
concat.write_phylip_to_file(sys.argv[3])


