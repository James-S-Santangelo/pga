# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies

rule quast_psic_ref:
    input:
        fasta = PSIC_REFERENCE, 
    output:
        directory(f"{QC_DIR}/quast/psic_ref")
    log: LOG_DIR + '/quast/quast_psic_ref.log'
    conda: '../envs/qc.yaml'
    threads: 8
    shell:
        """
        quast.py -o {output} \
            --threads {threads} \
            --split-scaffolds \
            --large \
            --k-mer-stats \
            {input.fasta} &> {log}
        """

# rule run_busco:
#     input:
#         prot = f"{rules.funannotate_annotate.output}/annotate_results/Trifolium_repens.proteins.fa"
#     output:
#         directory(f"{QC_DIR}/busco/TrR_v6_{{db}}")
#     log: LOG_DIR + '/busco/busco_{db}.log'
#     conda: '../envs/qc.yaml'
#     threads: 32
#     params:
#         out_path = f"{QC_DIR}/busco/",
#         out_name = "TrR_v6_{db}"
#     shell:
#         """
#         busco -m protein \
#             -i {input.prot} \
#             -o {params.out_name} \
#             --out_path {params.out_path} \
#             --lineage {wildcards.db} \
#             --force \
#             --cpu {threads} &> {log}
#         """

rule qc_done:
    input:
        rules.quast_haploid_ref.output,
        #expand(rules.run_busco.output, db = ['embryophyta_odb10', 'fabales_odb10'])
    output:
        f'{QC_DIR}/qc.done'
    shell:
        """
        touch {output}
        """
