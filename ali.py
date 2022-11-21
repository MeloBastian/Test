#! /usr/bin/python3
import sys
import os
import re
import random

univ_trans = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_codon(codon):
    if codon in univ_trans:
        return univ_trans[codon]
    else:
        return "?"

def translate_nucseq(nucseq):
    nsite = len(nucseq)
    if nsite % 3:
        print("error in translate: not multiple of 3")
        sys.exit()
    naa = nsite // 3
    return "".join([translate_codon(nucseq[3*i:3*(i+1)]) for i in range(naa)])

def get_gc3(nucseq):
    nsite = len(nucseq)
    if nsite % 3:
        print("error in translate: not multiple of 3")
        sys.exit()
    naa = nsite // 3
    ngc = 0
    ntot = 0
    for i in range(naa):
        if nucseq[3*i+2] in "gcGC":
            ngc += 1
        if nucseq[3*i+2] in "atgcATGC":
            ntot += 1
    f = 0
    if ntot:
        f = ngc / ntot
    return f

class SequenceAlignment:

    def __init__(self, file_name = "", file_stream = "", ali = "", format = "phylip", fromdict = dict()):
        self.nsite = 0
        self.ali = dict()
        if file_name:
            if format == "phylip":
                self.read_sequential_phylip_from_file(file_name)
            if format == "fasta":
                self.read_fasta_from_file(file_name)

        if file_stream:
            if format != "phylip":
                print("error in SequenceAlignment: stream constructor only accepts phylip format")
                sys.exit()
            self.read_sequential_phylip_from_stream(file_stream)

        if len(fromdict):
            self.ali = fromdict
            self.nsite = max([len(seq) for tax,seq in self.ali.items()])

    def get_taxon_set(self):
        return {tax for (tax,seq) in self.ali.items()}

    def get_ntaxa(self):
        return len(self.ali)

    def get_total_size(self):
        return sum([len(seq) for (tax,seq) in self.ali.items()])

    def get_nsite(self):
        return self.nsite

    def get_seq(self, tax):
        if tax not in self.ali:
            print("error in SequenceAlignment.get_seq: taxon not in alignment")
            sys.exit(1)
        return self.ali[tax]

    def read_sequential_phylip_from_file(self, filename):
        with open(filename, 'r') as infile:
            self.read_sequential_phylip_from_stream(infile)

    def read_sequential_phylip_from_stream(self, infile):
        (ntax,npos) = infile.readline().rstrip('\n').split()
        ntaxa = int(ntax)
        self.nsite = int(npos)

        self.ali = dict()
        for i in range(ntaxa):
            line = infile.readline()
            if not line:
                print("error when reading phylip")
                sys.exit()
            pair = line.rstrip('\n').split()
            if len(pair) != 2:
                print("error in readphylip: should have 2 fields per line")
                print(pair)
                sys.exit()
            if len(pair[1]) != self.nsite:
                print("error in read_phylip: non matching number of sites")
                sys.exit()
            self.ali[pair[0]] = pair[1]

    def read_fasta_from_file(self, filename):
        with open(filename, 'r') as infile:
            self.read_fasta_from_stream(infile)

    def read_fasta_from_stream(self, infile):
        taxname = ''
        for line in infile:
            line = line.rstrip('\n')
            if line:
                if line[0] == '>':
                    taxname = line[1:len(line)]
                    self.ali[taxname] = ''
                else:
                    self.ali[taxname] = self.ali[taxname] + line
        self.nsite = 0
        for (tax,seq) in self.ali.items():
            if not self.nsite:
                self.nsite = len(seq)
            else:
                pass
            #    if self.nsite != len(seq):
            #        print("error in read_fasta: sequences do not have same length")
            #        sys.exit()

    def prune_all_missing_taxa(self, noninformative_characters = "?X-*"):
        taxlist = []
        for (tax,seq) in self.ali.items():
            informative  = sum([s not in noninformative_characters for s in seq])
            if not informative:
                taxlist.append(tax)
        for tax in taxlist:
            del self.ali[tax]

    def countmin(self, min):
        n = 0
        for (tax,seq) in self.ali.items():
            if len(seq) <= min:
                n = n + 1
        return n
        
    def prune_min_length(self, min):
        newali = dict()
        for (tax,seq) in self.ali.items():
            if len(seq) > min:
                newali[tax] = seq
        del self.ali
        self.ali = newali

    def prune_gappy_taxa(self, minfrac, noninf = "?X-*"):
        cutoff = int(minfrac * self.get_nsite())
        taxlist = []
        for (tax,seq) in self.ali.items():
            if len([s for s in seq if s not in noninf]) < min:
                taxlist.append(tax)
        for tax in taxlist:
            del self.ali[tax]

    def prune_gappy_columns(self, minfrac, noninf = "?X-*"):
        nmissing = [sum([seq[i] in noninf for tax,seq in self.ali.items()]) for i in range(self.get_nsite())]
        cutoff = int(minfrac * self.get_ntaxa())
        new_nsite = sum([m <= cutoff for m in nmissing])
        for (tax,seq) in self.ali.items():
            newseq = "".join([seq[i] for i in range(len(nmissing)) if nmissing[i] <= cutoff])
            self.ali[tax] = newseq
        self.nsite = new_nsite

    def prune_all_gappy_columns(self, noninf = "?X-*"):
        nmissing = [sum([seq[i] in noninf for tax,seq in self.ali.items()]) for i in range(self.get_nsite())]
        new_nsite = sum([m < self.get_ntaxa() for m in nmissing])
        for (tax,seq) in self.ali.items():
            newseq = "".join([seq[i] for i in range(len(nmissing)) if nmissing[i] < self.get_ntaxa()])
            self.ali[tax] = newseq
        self.nsite = new_nsite

    def unalign(self):
        for (tax,seq) in self.ali.items():
            newseq = "".join([s for s in seq if s not in "-"])
            self.ali[tax] = newseq

    def redefine_taxon_set(self, taxon_set, all_taxa = True, noninformative_characters = "?X-*"):
        temp_tax_list = [tax for (tax,seq) in self.ali.items()]
        for tax in temp_tax_list:
            if tax not in taxon_set:
                del self.ali[tax]
        if all_taxa:
            all_missing  = "?" * self.nsite
            for tax in taxon_set:
                if tax not in self.ali:
                    self.ali[tax] = all_missing
        else:
            self.prune_all_missing_taxa(noninformative_characters = noninformative_characters)

    def taxon2genus(self):
        match_genus = r"^([A-Z][a-z]+)\_.*$"
        temp_tax_list = [tax for (tax,seq) in self.ali.items()]
        newali = dict()
        for tax in temp_tax_list:
            m = re.match(match_genus, tax)
            if not m:
                print("error in taxon2genus: could not parse taxon name {0}".format(tax))
                sys.exit(1)
            genus = m.group(1)
            if genus in newali:
                print("in taxon2genus: genus {0} appears multiple times".format(genus))
                print("over-writing previous sequence")
            newali[genus] = self.ali[tax]

        for tax in temp_tax_list:
            del self.ali[tax]

        for (genus,seq) in newali.items():
            self.ali[genus] = seq

    def write_phylip_to_file(self, filename, force = False):
        if not force and os.path.exists(filename):
            print("in write_phylip: destination file already exists: ", filename)
            sys.exit()
        with open(filename, 'w') as outfile:
            self.write_phylip_to_stream(outfile)

    def write_randomized_phylip_to_file(self, filename, force = False):
        if not force and os.path.exists(filename):
            print("in write_phylip: destination file already exists: ", filename)
            sys.exit()
        with open(filename, 'w') as outfile:
            self.write_randomized_phylip_to_stream(outfile)

    def write_phylip_to_stream(self, outfile):
        outfile.write("{0} {1}\n".format(self.get_ntaxa(), self.nsite))
        for tax in sorted(self.ali.keys()):
        # for (tax,seq) in self.ali.items():
            seq = self.ali[tax]
            outfile.write("{0}  {1}\n".format(tax,seq))

    def write_randomized_phylip_to_stream(self, outfile):
        outfile.write("{0} {1}\n".format(self.get_ntaxa(), self.nsite))
        siteranks = random.sample([i for i in range(self.get_nsite())], self.get_nsite())
        for tax in sorted(self.ali.keys()):
            seq = self.ali[tax]
            outfile.write("{0}  {1}\n".format(tax,"".join([seq[i] for i in siteranks])))

    def write_nexus_to_file(self, filename, force = False):
        if not force and os.path.exists(filename):
            print("in write_phylip: destination file already exists: ", filename)
            sys.exit()
        with open(filename, 'w') as outfile:
            self.write_nexus_to_stream(outfile)

    def write_nexus_to_stream(self, outfile):
        outfile.write("#NEXUS\n")
        outfile.write("\n")
        outfile.write("BEGIN DATA;\n")
        outfile.write("DIMENSIONS NTAX={0} NCHAR={1};\n".format(self.get_ntaxa(), self.get_nsite()))
        outfile.write("FORMAT DATATYPE=DNA GAP=- MISSING=?;\n")
        outfile.write("MATRIX\n")
        outfile.write("\n")
        for (tax,seq) in self.ali.items():
            outfile.write("{0}  {1}\n".format(tax,seq))
        outfile.write("\t;\n")
        outfile.write("end;\n")

    def write_fasta_to_file(self, filename, force = False):
        if not force and os.path.exists(filename):
            print("in write_fasta: destination file already exists: ", filename)
            sys.exit()
        with open(filename, 'w') as outfile:
            self.write_fasta_to_stream(outfile)

    def write_fasta_to_stream(self, outfile):
        for (tax,seq) in self.ali.items():
            outfile.write(">{0}\n{1}\n".format(tax,seq))

    def sub_alignment(self, mask):
        ali = dict()
        for (tax, seq) in self.ali.items():
            newseq = ''.join([c for (i,c) in enumerate(seq) if mask[i]])
            ali[tax] = newseq
        newali = SequenceAlignment()
        newali.nsite = sum(mask)
        newali.ali = ali
        return newali

    def split(self, mask):
        ali1 = dict()
        ali2 = dict()
        for (tax, seq) in self.ali.items():
            newseq1 = ''.join([c for (i,c) in enumerate(seq) if mask[i]])
            newseq2 = ''.join([c for (i,c) in enumerate(seq) if not mask[i]])
            ali1[tax] = newseq1
            ali2[tax] = newseq2
        newali1 = SequenceAlignment()
        newali1.nsite = sum(mask)
        newali1.ali = ali1
        newali2 = SequenceAlignment()
        newali2.nsite = self.nsite - sum(mask)
        newali2.ali = ali2
        return (newali1, newali2)

    def split_equal(self, n):
        size = int(self.nsite / n)
        print(size)
        alis = []
        for k in range(n):
            ali1 = dict()
            for (tax, seq) in self.ali.items():
                newseq = seq[k*size:(k+1)*size]
                ali1[tax] = newseq
            newali = SequenceAlignment()
            newali.nsite = size
            newali.ali = ali1
            alis.append(newali)
        return alis

    def translate(self, code = univ_trans):
        if self.nsite % 3:
            print("error in translate: not multiple of 3")
            sys.exit()
        naa = self.nsite // 3
        ali = dict()
        for (tax,nucseq) in self.ali.items():
            ali[tax] = translate_nucseq(nucseq)
        protali = SequenceAlignment()
        protali.nsite = naa
        protali.ali = ali
        return protali

    def get_mean_diversity(self, noninformative_characters = "?X-*"):
        sets = [set() for i in range(self.nsite)]
        for tax,seq in self.ali.items():
            for i in range(self.nsite):
                if seq[i] not in noninformative_characters:
                    sets[i].add(seq[i])
        return sum([len(s) for s in sets])/len(sets)

    def get_gc3(self):
        meangc3 = 0
        for tax,seq in self.ali.items():
            meangc3 += get_gc3(seq)
        meangc3 /= len(self.ali)
        return meangc3

class MultiGeneSequenceAlignment:

    def __init__(self, dir_name = "", list_name = "", file_name = "", alignments = "", format = "phylip", header=False):

        self.alignments = dict()
        self.taxon_set = set()

        if list_name:
            self.read_from_list(dir_name, list_name, format = format, header= header)
        if file_name:
            self.read_from_file(file_name, format = format)
        if alignments:
            self.alignments = alignments

        self.make_taxon_set()

    def make_taxon_set(self):
        self.taxon_set = set()
        for (gene,ali) in self.alignments.items():
            self.taxon_set = self.taxon_set | ali.get_taxon_set()

    def get_taxon_set(self):
        return self.taxon_set

    def get_ntaxa(self):
        return len(self.taxon_set)

    def prune_min_ntaxa(self, min):
        genelist = []
        for (gene,ali) in self.alignments.items():
            if ali.get_ntaxa() < min:
                genelist.append(gene)
        for gene in genelist:
            del self.alignments[gene]

    def prune_min_nsite(self, min):
        genelist = []
        for (gene,ali) in self.alignments.items():
            if ali.get_nsite() < min:
                genelist.append(gene)
        for gene in genelist:
            del self.alignments[gene]

    def prune_gappy_columns(self, minfrac, noninf = "?X-*"):
        for gene,ali in self.alignments.items():
            ali.prune_gappy_columns(minfrac, noninf=noninf)

    def prune_all_gappy_columns(self, noninf = "?X-*"):
        for gene,ali in self.alignments.items():
            ali.prune_all_gappy_columns(noninf=noninf)

    def concatenate(self):
        ali = dict()
        for taxon in self.taxon_set:
            seq = "".join([self.alignments[gene].get_seq(taxon) for gene in self.alignments])
            ali[taxon] = seq
        return SequenceAlignment(fromdict=ali)

    def homogeneize_taxon_sets(self):
        for (gene,ali) in self.alignments.items():
            ali.redefine_taxon_set(self.taxon_set, all_taxa=True)

    def get_ngene(self):
        return len(self.alignments)

    def read_from_list(self, dir_name, list_name, format = "phylip", header = False):
        with open(list_name, 'r') as infile:
            ngene = 0
            if header:
                ngene = int(infile.readline().rstrip('\n'))
            gene_list = [line.rstrip('\n') for i,line in enumerate(infile) if not ngene or i<ngene]
            for gene in gene_list:
                ali = SequenceAlignment(dir_name + "/" + gene, format = format)
                self.alignments[gene] = ali

    def read_from_file(self, file_name, format = "phylip"):
        with open(file_name, 'r') as infile:
            header = infile.readline().rstrip('\n')
            if header != "ALI":
                print("error in MultiGeneSequenceAlignment.read_from_file: header")
                sys.exit()

            ngene = int(infile.readline().rstrip('\n'))

            for i in range(ngene):
                gene = infile.readline().rstrip('\n')
                ali = SequenceAlignment()
                ali.read_sequential_phylip_from_stream(infile)
                self.alignments[gene] = ali
                
    def write_all_to_file(self, filename):
        with open(filename, 'w') as outfile:
            outfile.write("ALI\n")
            outfile.write("{0}\n".format(self.get_ngene()))
            for (gene,ali) in self.alignments.items():
                outfile.write("{0}\n".format(gene))
                ali.write_phylip_to_stream(outfile)

    def write_all_to_files(self, basename, force = False):
        for gene,ali in self.alignments.items():
            filename = basename + gene
            if not force and os.path.exists(filename):
                print("error when writing gene-specific alignments: files {0} already exists".format(filename))
                sys.exit(1)
            with open(filename, 'w') as outfile:
                ali.write_phylip_to_stream(outfile)

    def write_genes_to_nexus(self, dirname):
        for (gene,ali) in self.alignments.items():
            ali.write_nexus_to_file(dirname + "/" + gene)

    def get_gene_list(self):
        gene_list = sorted([gene for (gene,ali) in self.alignments.items()])
        return gene_list

    def get_gene_ali(self, gene):
        if gene not in self.alignments:
            print("error: gene not in multigene sequence alignment")
            sys.exit(1)
        return self.alignments[gene]

    def change_gene_names(self, name_pattern):
        alignments = dict()
        for (gene, ali) in self.alignments:
            m = re.match(pattern, gene)
            if not m:
                print("error in MultiGeneSequenceAlignment.read_from_list: gene name does not match pattern")
                sys.exit()
            new_name = m.group(1)
            alignments[new_name] = self.alignments[gene]
        self.alignments = alignments

    def prune_all_missing_taxa(self, noninformative_characters = "?X-*"):
        for (gene,ali) in self.alignments.items():
            ali.prune_all_missing_taxa(noninformative_characters = noninformative_characters)

    def redefine_taxon_set(self, taxon_set, all_taxa = True, noninformative_characters = "?X-*"):
        for (gene,ali) in self.alignments.items():
            ali.redefine_taxon_set(taxon_set, all_taxa = all_taxa, noninformative_characters = noninformative_characters)
        self.make_taxon_set()

    def translate(self, code = univ_trans):
        alignments = dict()
        for (gene,ali) in self.alignments.items():
            protali = ali.translate(code = code)
            alignments[gene] = protali
        prot_multiali = MultiGeneSequenceAlignment(alignments=alignments)
        return prot_multiali

    def get_mean_diversity(self, gene, noninformative_characters="?X-*"):
        return self.alignments[gene].get_mean_diversity(noninformative_characters = noninformative_characters)

    def get_gc3(self, gene):
        return self.alignments[gene].get_gc3()

