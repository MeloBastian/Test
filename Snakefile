#   snakemake --cluster "sbatch -J {params.name} -p normal -N 1 --ntasks={params.cpus} --mem={params.mem} --constraint='skylake|haswell|broadwell'  -t {params.time} -o {params.out} -e {params.err}" -j 30 --rerun-incomplete -n


#list of sp and genes selected
PathData="/beegfs/data/mbastian/Zoonomia_postbusco/"

#sp_list_path=PathData+"/Genes_to_analyse/180splist"
sp_list_path=PathData+"/Genes_to_analyse_V2/180sp"
sp_file = open(sp_list_path, "r")
#splist = [line.rstrip("\n") for line in spfile]

#gene_list_path=PathData+"/Genes_to_analyse/7544genelist"
gene_list_path=PathData+"/Genes/Genes_V2/genelist_filtred_V2"
gene_file=open(gene_list_path,"r")
genelist=gene_file.readlines()
genelist = [line.rstrip("\n") for line in genelist]
#genelist=genelist[0:20]


localrules: filter_sp, conc_and_align_1genus, conc_and_ali_allsp
#print(len(genelist))


rule all:
    input:
        expand(PathData +"Genes_to_analyse_V2/{gene}_filtred_ali.fna", gene=genelist), #prank
        expand(PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_aa.fna", gene=genelist), #seaview
        expand(PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_aa_hmm.fasta", gene=genelist),#hmm
        expand(PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_aa_hmm.doneflag",gene=genelist),#hmm
        expand(PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_NT.doneflag",gene=genelist),  #macse
        expand(PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_hmm_bmge.fna", gene=genelist), #bmge
        expand(PathData+ "Genes_to_analyse_V2/{gene}_filtred_ali_hmm_bmge.doneflag", gene=genelist) #bmge

        #expand( PathData + "Genes_to_analyse/Zoonomia_Busco_180sp.conc")

rule filter_sp:
    output:
        expand(PathData + "Genes_to_analyse_V2/{gene}_filtred.fna",gene=genelist)
    shell:
        """
        python -u removespfromali.py {gene_list_path}
        """

rule Prank:
    params:
        time="40:00:00", mem="3G", name="Prank", cpus="1",
        out=PathData + "Log_Out/Prank/{gene}_prank_out",
        err=PathData + "Log_Out/Prank/{gene}_prank_err"
    input :
        file = PathData + "Genes_to_analyse_V2/{gene}_filtred.fna"
    output:
        gene_pranked = PathData +"Genes_to_analyse_V2/{gene}_filtred_ali.fna"
    shell:
        """
        echo $SLURM_JOB_NODELIST
        #conda activate /beegfs/home/mbastian/bioinfo
        prank -d={input.file} -o={output.gene_pranked} -DNA -codon
        mv {output.gene_pranked}.best.fas {output.gene_pranked}
        """

rule translate:
    params:
        time="00:10:00", mem="3G", name="Translation", cpus="1",
        out=PathData + "Log_Out/seaview_tr/{gene}_seaview_tr.out",
        err=PathData + "Log_Out/seaview_tr/{gene}_seaview_tr.err",constraint="skylake|haswell|broadwell",
        #dir= PathData+"Genes_to_analyse/"
    input:
        file= PathData+"Genes_to_analyse_V2/{gene}_filtred_ali.fna"
    output:
        genes_tr = PathData+ "Genes_to_analyse_V2/{gene}_filtred_ali_aa.fna"
    shell:
        """
        /beegfs/data/anguyentr/tools/seaview/seaview -convert -translate -o {output.genes_tr} {input.file}
        """


rule Hmm_cleaner:
    params:
        time="1:00:00", mem="3G", name="HmmCl", cpus="1",
        out=PathData + "Log_Out/HmmCl/{gene}_hmm.out",
        err=PathData + "Log_Out/HmmCl/{gene}_hmm.err",constraint="skylake|haswell|broadwell",
        file= "{gene}_filtred_ali_aa.fna",
        dir = PathData+"Genes_to_analyse_V2/" #hmmCl can't read a path
    input:
        file=  PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_aa.fna"
    output:
        genes_hmm = PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_aa_hmm.fasta",
        #out=PathData+"Genes_to_analyse/{gene}_filtred_ali_hmm.fasta",
        done= PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_aa_hmm.doneflag"
    shell:
        """
        cd {params.dir}
        echo $SLURM_JOB_NODELIST
        singularity exec /beegfs/home/tricou/singularity_image/hmm.sif HmmCleaner.pl {params.file}
        
        echo done > {output.done}
        """

rule Macse:
    params:
        time="00:10:00", mem="3G", name="Macse", cpus="1",
        out=PathData + "Log_Out/Macse/{gene}_macse.out",
        err=PathData + "Log_Out/Macse/{gene}_macse.err",constraint="skylake|haswell|broadwell"
    input:
        fileaa= PathData+"Genes_to_analyse_V2/{gene}_filtred_ali_aa_hmm.fasta",
        filenuc= PathData +"Genes_to_analyse_V2/{gene}_filtred_ali.fna"
    output:
        genes_macse = PathData +"Genes_to_analyse_V2/{gene}_filtred_ali_NT.fna",
        done= PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_NT.doneflag"
    shell:
        """
        /beegfs/home/mbastian/bioinfo/bin/java -jar /beegfs/data/mbastian/macse_v2.06.jar -prog reportMaskAA2NT -align_AA {input.fileaa}  -align {input.filenuc} -mask_AA -
        echo done > {output.done}
        """

rule BMGE:
    params:
        time="01:00:00", mem="3G", name="BMGE", cpus="1",
        out=PathData + "Log_Out/BMGE/{gene}_prank_hmm_bmge_out",
        err=PathData + "Log_Out/BMGE/{gene}_prank_hmm_bmge_err",constraint="skylake|haswell|broadwell",
        dir= PathData+"Genes_to_analyse_V2/"
    input :
        file =  PathData +"Genes_to_analyse_V2/{gene}_filtred_ali_NT.fna"
    output:
        gene_bmge = PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_hmm_bmge.fna",
        done= PathData + "Genes_to_analyse_V2/{gene}_filtred_ali_hmm_bmge.doneflag"
    shell:
        """
        /beegfs/home/mbastian/bioinfo/bin/java -jar /beegfs/data/mbastian/BMGE.jar -i {input.file} -t CODON -o {output.gene_bmge}
        echo done > {output.done}
        """

    # rule conc_and_align_1genus: #create a conc and multiali for the sp in splist_1genus
    #     input:
    #         file=PathData + "Genes_to_analyse/{gene}_filtred_ali_hmmCl_bmge.fna"
    #     output:
    #         gene_multiali= PathData + "Genes_to_analyse/Zoonomia_Busco_Xsp.ali",
    #         gene_conc= PathData + "Genes_to_analyse/Zoonomia_Busco_Xsp.conc"
    #     shell:
    #         """
    #         python -u ali_and_conc_zoo_1genus.py genelist {output.gene_multiali} {output.gene_conc}
    #         """

    # rule conc_and_align_allsp:
    #     input:
    #         file=PathData + "Genes_to_analyse/{gene}_filtred_ali_hmmCl_bmge.fna"
    #     output:
    #         gene_multiali= PathData + "Genes_to_analyse/Zoonomia_Busco_180sp.ali",
    #         gene_conc= PathData + "Genes_to_analyse/Zoonomia_Busco_180sp.conc"
    #     shell:
    #         """
    #         python -u ali_and_conc_zoo_180sp.py genelist {output.gene_multiali} {output.gene_conc}
    #         """
