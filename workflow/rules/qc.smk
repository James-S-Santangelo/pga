# Rule to QC original Dovetail haplotypes, revised haplotypes, and final assemblies

rule quast_psic_ref:
    """
    Run QUAST to generate genome assembly statistics
    """
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
    """
    Run BUSCO in genome-mode against tetrapod database
    """
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

rule functional_stats:
    """
    Generate summary statistics of functional annotation using AGAT
    """
    input:
        gff = rules.gff_sort_functional.output,
        fasta = rules.repeat_masker.output.fasta
    output:
        directory(f"{QC_DIR}/agat_funcStats")
    log: f"{LOG_DIR}/agat_funcStats/agat_funcStats.log"
    container: 'docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0'
    shell:
        """
        agat_sp_functional_statistics.pl --gff {input.gff} \
                -gs {input.fasta} \
                --output {output} 2> {log}
        """

rule run_busco_protein:
    """
    Run BUSCO in protein-mode against tetrapod database
    """
    input:
        rules.get_proteins_finalGFF.output
    output:
        directory(f"{QC_DIR}/busco/protein/psic_protein_tetrapoda")
    log: LOG_DIR + '/busco/busco_protein_tetrapoda.log'
    conda: '../envs/qc.yaml'
    threads: 32
    params:
        out_path = f"{QC_DIR}/busco/protein/",
        out_name = "psic_protein_tetrapoda"
    shell:
        """
        busco -m protein \
                -i {input} \
                -o {params.out_name} \
                --out_path {params.out_path} \
                --lineage tetrapoda_odb10 \
                --force \
                --cpu {threads} &> {log}
        """


rule qc_done:
    input:
        rules.quast_psic_ref.output,
        rules.run_busco_genome.output,
        rules.functional_stats.output,
        rules.run_busco_protein.output
    output:
        f'{QC_DIR}/qc.done'
    shell:
        """
        touch {output}
        """
