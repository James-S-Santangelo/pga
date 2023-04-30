# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies

rule quast_psic_ref:
    input:
        fasta = rules.setup_ref.output.fasta 
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

rule run_busco_genome:
    input:
        fasta = rules.subset_ref_bySeqLength.output.fasta
    output:
        directory(f"{QC_DIR}/busco/genome/psic_genome_tetrapoda")
    log: LOG_DIR + '/busco/busco_genome_tetrapoda.log'
    conda: '../envs/qc.yaml'
    threads: 32
    params:
        out_path = f"{QC_DIR}/busco/genome/",
        out_name = "psic_genome_tetrapoda"
    shell:
        """
        busco -m genome \
            -i {input.fasta} \
            -o {params.out_name} \
            --out_path {params.out_path} \
            --lineage tetrapoda_odb10 \
            --force \
            --cpu {threads} &> {log}
        """

rule qc_done:
    input:
        rules.quast_psic_ref.output,
        rules.run_busco_genome.output
    output:
        f'{QC_DIR}/qc.done'
    shell:
        """
        touch {output}
        """
