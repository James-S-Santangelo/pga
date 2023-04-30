rule setup_ref:
    input:
        fasta = PSIC_REF
    output:
        fasta = f"{REF_DIR}/psic_genome/psic_refererence.fasta",
        fai = f"{REF_DIR}/psic_genome/psic_refererence.fasta.fai"
    conda: '../envs/ref.yaml'
    shell:
        """
        ln -sf {input} {output.fasta} &&
        samtools faidx {output.fasta}
        """


rule subset_ref_bySeqLength:
    input:
        fasta = rules.setup_ref.output.fasta,
        fai = rules.setup_ref.output.fai
    output:
        fasta = f"{REF_DIR}/psic_genome/psic_ref_100Kb.fasta",
        scaffs = f"{PROGRAM_RESOURCE_DIR}/scaffs_over100Kb.txt"
    conda: '../envs/ref.yaml'
    shell:
        """
        cat {input.fai} | awk '$2 > 100000 {{print $1}}' > {output.scaffs}&&
        samtools faidx {input.fasta} -r {output.scaffs} > {output.fasta}
        """
