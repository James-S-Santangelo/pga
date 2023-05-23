# Rules to annotate the genome

#####################################
#### REPEAT MODELING AND MASKING ####
#####################################

rule configure_repbase:
    """
    Configure RepBase Database for use with RepeatModeler
    """
    input:
        REPBASE
    output:
        libdir = directory(f"{PROGRAM_RESOURCE_DIR}/Libraries"),
        dfam = f"{PROGRAM_RESOURCE_DIR}/Libraries/Dfam.h5",
        rm = f"{PROGRAM_RESOURCE_DIR}/Libraries/RepeatMaskerLib.h5"
    container: 'docker://dfam/tetools:1.8'
    log: LOG_DIR + '/configure_repbase/configure_repbase.log'
    params:
        untar_dir = f"{PROGRAM_RESOURCE_DIR}"
    shell:
        """
        ( tar -xzf {input} -C {params.untar_dir} &&
         cp -r /opt/RepeatMasker/Libraries/* {output.libdir} &&
         ln -sf {output.dfam} {output.rm} &&
         addRepBase.pl --libdir {output.libdir} ) &> {log}
        """

rule build_repeat_modeler_db:
    """
    Build RepeatModeler Database from haploid mapping reference assembly
    """
    input:
        rules.subset_ref_bySeqLength.output.fasta
    output:
        multiext(f"{ANNOTATION_DIR}/repeat_modeler/rmdb", '.nhr', '.nin', '.nnd', '.nni', '.nog', '.nsq', '.translation')
    container: 'docker://dfam/tetools:1.8'
    log: LOG_DIR + '/build_repeat_modeler_db/rmdb.log'
    params:
        out = f"{ANNOTATION_DIR}/repeat_modeler/rmdb"
    shell:
        """
        BuildDatabase -name {params.out} {input} &> {log}
        """

rule repeat_modeler:
    """
    Build de novo repeat library from haploid mapping reference assembly using RepeatModeler
    """
    input:
        rules.build_repeat_modeler_db.output,
        rules.configure_repbase.output
    output:
        fasta = f"{ANNOTATION_DIR}/repeat_modeler/rmdb-families.fa",
        stk = f"{ANNOTATION_DIR}/repeat_modeler/rmdb-families.stk"
    threads: 24
    container: 'docker://dfam/tetools:1.8'
    log: LOG_DIR + '/repeat_modeler/rm.log'
    params:
        db_base = lambda wildcards, input: os.path.splitext(input[0])[0]
    shell:
        """
        RepeatModeler -database {params.db_base} \
            -threads {threads} \
            -LTRStruct &> {log}
        """

rule merge_repeat_databases:
    """
    Merge de novo repeat library with Green Plant (i.e., Viridiplantae) repeat library from RepBase
    """
    input:
        rm_db = rules.configure_repbase.output.rm,
        tr_db = rules.repeat_modeler.output.fasta
    output:
        f"{PROGRAM_RESOURCE_DIR}/Libraries/rm_merged_db.fasta"
    container: 'docker://dfam/tetools:1.8'
    log: LOG_DIR + '/merge_repeat_databases/merge_repeat_databases.log'
    params:
        rm_db_fasta = f"{PROGRAM_RESOURCE_DIR}/Libraries/rm_db.fasta"
    shell:
        """
        ( famdb.py -i {input.rm_db} families \
            --format fasta_name \
            --ancestors --descendants --include-class-in-name \
            'tetrapoda' > {params.rm_db_fasta} &&

        cat {params.rm_db_fasta} {input.tr_db} > {output} ) 2> {log}
        """

rule repeat_masker:
    """
    Softmask repeats in the white clover haploid mapping reference assembly.
    """
    input:
        lib = rules.merge_repeat_databases.output,
        fasta = rules.subset_ref_bySeqLength.output.fasta
    output:
        fasta = f"{ANNOTATION_DIR}/repeat_masker/psic_reference_softMasked.fasta",
        cat = f"{ANNOTATION_DIR}/repeat_masker/psic_reference_repeatMasker.cat.gz",
        out = f"{ANNOTATION_DIR}/repeat_masker/psic_reference_repeatMasker.out",
        gff = f"{ANNOTATION_DIR}/repeat_masker/psic_reference_repeatMasker.gff",
        stats = f"{ANNOTATION_DIR}/repeat_masker/psic_reference_repeatMasker.tbl"
    threads: 24
    container: 'docker://dfam/tetools:1.8'
    log: LOG_DIR + '/repeat_masker/repeat_masker.log'
    params:
        outdir = f"{ANNOTATION_DIR}/repeat_masker/",
        gc = 44
    shell:
        """
        ( RepeatMasker -e ncbi \
            -pa {threads} \
            -xsmall \
            -gc {params.gc} \
            -lib {input.lib} \
            -dir {params.outdir} \
            -gff {input.fasta} &&
            
            mv {params.outdir}/*.fasta.masked {output.fasta}
            mv {params.outdir}/*.cat.gz {output.cat}
            mv {params.outdir}/*.out {output.out}
            mv {params.outdir}/*.gff {output.gff}
            mv {params.outdir}/*.tbl {output.stats} ) &> {log}
        """

###########################################
#### DOWNLOAD AND MAP ALL RNASEQ READS ####
###########################################

rule prefetch:
    """
    Pre-fetch RNAseq libraries that will be used for structural gene annotation
    """
    output:
        temp(directory(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}"))
    log: LOG_DIR + '/prefetch/{acc}.log'
    conda: '../envs/annotation.yaml'
    params:
        outdir = f"{ANNOTATION_DIR}/rnaseq_reads"
    shell:
        """
        prefetch {wildcards.acc} -O {params.outdir} 2> {log}
        """

rule fasterq_dump:
    """
    Download pre-fetched RNAseq libraries used for structural gene annotation
    """
    input:
        expand(rules.prefetch.output, acc=RNASEQ_ACCESSIONS)
    output:
        R1 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_1.fastq"),
        R2 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_2.fastq")
    log: LOG_DIR + '/fastq_dump/{acc}_fastq_dump.log'
    conda: '../envs/annotation.yaml'
    params:
        outdir = f"{ANNOTATION_DIR}/rnaseq_reads"
    threads: 6
    shell:
        """
        cd {params.outdir}
        fasterq-dump --split-3 \
            -e {threads} \
            --skip-technical \
            {wildcards.acc} 2> {log}
        """

rule gzip_fastq:
    """
    Gzip downloaded RNAseq libraries
    """
    input:
        rules.fasterq_dump.output
    output:
        R1 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_1.fastq.gz"),
        R2 = temp(f"{ANNOTATION_DIR}/rnaseq_reads/{{acc}}_2.fastq.gz")
    log: LOG_DIR + '/gzip_fastq/{acc}_gzip.log'
    shell:
        """
        gzip {input} 2> {log} 
        """

rule build_star:
    """
    Build STAR Database from softmasked haploid mapping reference assembly
    """
    input:
        masked_genome = rules.repeat_masker.output.fasta 
    output:
        temp(directory(f"{ANNOTATION_DIR}/star/star_build"))
    log: LOG_DIR + '/star/star_build.log'
    conda:'../envs/annotation.yaml'
    threads: 6
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.masked_genome} \
            --runThreadN {threads} &> {log}
        """

rule align_star:
    """
    Align downloaded RNAseq reads to softmasked haploid reference assembly using STAR in two-pass mode.
    Set max intron length to 10000.
    """
    input:
        star_build = rules.build_star.output,
        R1 = rules.gzip_fastq.output.R1,
        R2 = rules.gzip_fastq.output.R2
    output:
        star_align = temp(f"{ANNOTATION_DIR}/star/star_align/{{acc}}/{{acc}}_Aligned.sortedByCoord.out.bam") 
    log: LOG_DIR + '/star/{acc}_star_align.log'
    conda:'../envs/annotation.yaml'
    params:
        out = f"{ANNOTATION_DIR}/star/star_align/{{acc}}/{{acc}}_"
    threads: 6
    shell:
        """
        STAR --readFilesIn {input.R1} {input.R2} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --twopassMode Basic \
            --genomeDir {input.star_build} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --alignIntronMax 10000 \
            --outFileNamePrefix {params.out} &> {log} 
        """

rule merge_rnaseq_bams:
    """
    Merges STAR-aligned RNAseq reads into a single BAM file. This will be used as input to BRAKER in RNAseq-mode
    """
    input:
        Star_Bam = expand(rules.align_star.output, acc=RNASEQ_ACCESSIONS),
    output:
        bam = f"{ANNOTATION_DIR}/star/allBams_merged.bam",
        bai = f"{ANNOTATION_DIR}/star/allBams_merged.bam.bai"
    conda: '../envs/annotation.yaml'
    threads: 30
    log: LOG_DIR + '/merge_rnaseq_bams/bams_merge.log'
    shell:
        """
        ( samtools merge -@ {threads} -r -o {output.bam} {input} &&\
            samtools index {output.bam} ) 2> {log}
        """

###############################
###STRUCTURAL ANNOTATION ####
###############################

rule vertebrata_orthodb:
    """
    Download Vertebrate OrthoDB
    """
    output:
        ProrteinDB = f"{PROGRAM_RESOURCE_DIR}/orthodb/Vertebrata.fa"
    log: LOG_DIR + "/orthodb/ortho.log"
    params:
        outdir = f"{PROGRAM_RESOURCE_DIR}/orthodb"
    shell:
        """
        ( wget --no-check-certificate https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Vertebrata.fa.gz -P {params.outdir} && 
        gunzip {params.outdir}/Vertebrata.fa.gz ) 2> {log}
        """

rule download_podacris_proteins:
    """
    Download proteins from closely-related Podacris species (i.e., P. muralis and P. raffonei)
    """
    output:
        prot = f"{PROGRAM_RESOURCE_DIR}/podacris_proteins/{{psp}}/{{psp}}_proteins.fa"
    log: f"{LOG_DIR}/download_podacris_proteins/{{psp}}.log"
    params:
        ref_seq = lambda wildcards: 'GCF_027172205.1' if wildcards.psp == 'Praf' else 'GCF_004329235.1',
        outdir = f"{PROGRAM_RESOURCE_DIR}/podacris_proteins/{{psp}}/" 
    shell:
        """
        curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{params.ref_seq}/download?include_annotation_type=PROT_FASTA&filename={wildcards.psp}.zip" -H "Accept: application/zip" &&
        unzip {wildcards.psp}.zip -d {params.outdir} &&
        cp  {params.outdir}/ncbi_dataset/data/{params.ref_seq}/protein.faa {output}
        rm -rf {params.outdir}/ncbi_dataset
        """

rule download_uniprot_lacertidae_db:
    """
    Download all lacertidae proteins from UniProtKB
    """
    output:
        f'{PROGRAM_RESOURCE_DIR}/uniprot/lacertidae_proteins.fasta'
    log: LOG_DIR + '/uniprot_download/lacertidae_proteins_dl.log'
    params:
        url = "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28taxonomy_name%3Alacertidae%29"
    shell:
        """
        curl --output {output} '{params.url}' 2> {log}
        """

rule combine_protein_dbs:
    """
    Combine the vertebrate OrthoDB, Podacris, and lacertidae protein databases. This will be used as input to BRAKER in protein mode. 
    """
    input:
        ortho = rules.vertebrata_orthodb.output,
        psp = expand(rules.download_podacris_proteins.output, psp = ['Praf', 'Pmur']),
        unip = rules.download_uniprot_lacertidae_db.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/allProteins.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule braker_protein:
    """
    Run BRAKER in protein-mode
    """
    input:
        proteins = rules.combine_protein_dbs.output,
        masked_genome = rules.repeat_masker.output.fasta,
    output:
        prot_aug = f"{ANNOTATION_DIR}/braker/protein/Augustus/augustus.hints.gtf",
        hints_prot = f"{ANNOTATION_DIR}/braker/protein/hintsfile.gff"
    log: LOG_DIR + '/braker/proteins.log'
    params:
        outputdir = f"{ANNOTATION_DIR}/braker/protein"
    threads: 32
    container: 'docker://teambraker/braker3:v.1.0.4' 
    shell:
        """
        braker.pl --genome {input.masked_genome} \
            --prot_seq {input.proteins} \
            --softmasking \
            --useexisting \
            --threads {threads} \
            --workingdir {params.outputdir} \
            --species "Podacris siculus prot" 2> {log}
        """ 

rule braker_rnaseq:
    """
    Run BRAKER in RNAseq-mode
    """
    input:
        masked_genome = rules.repeat_masker.output.fasta,
        bam = rules.merge_rnaseq_bams.output.bam
    output:
        rna_aug = f"{ANNOTATION_DIR}/braker/rnaseq/Augustus/augustus.hints.gtf",
        hints_rna = f"{ANNOTATION_DIR}/braker/rnaseq/hintsfile.gff"
    log: LOG_DIR + '/braker/rnaseq.log'
    params:
        outputdir = f"{ANNOTATION_DIR}/braker/rnaseq"
    threads: 32
    container: 'docker://teambraker/braker3:v.1.0.4' 
    shell:
        """
        braker.pl --genome {input.masked_genome} \
            --bam {input.bam} \
            --softmasking \
            --useexisting \
            --threads {threads} \
            --workingdir {params.outputdir} \
            --species "Podacris siculus rna" 2> {log} 
        """

rule tsebra_combine:
    """
    Combine evidence from BRAKER Rnaseq and BRAKER protein using TSEBRA
    """
    input:
        rules.braker_protein.output,
        rules.braker_rnaseq.output,
        rna_aug = rules.braker_rnaseq.output.rna_aug, 
        prot_aug = rules.braker_protein.output.prot_aug, 
        hints_rna= rules.braker_rnaseq.output.hints_rna,
        hints_prot = rules.braker_protein.output.hints_prot 
    output:
        braker_combined = f"{ANNOTATION_DIR}/braker/tsebra/braker_combined.gtf"
    log: LOG_DIR + '/braker/tsebra.log'
    container: 'docker://teambraker/braker3:v.1.0.4'
    params:
        config = '../config/tsebra.cfg',
    shell:
        """
        tsebra.py -g {input.prot_aug} \
            -k {input.rna_aug} \
            -c {params.config} \
            -e {input.hints_rna},{input.hints_prot} \
            -o {output} 2> {log} 
        """ 

rule rename_tsebra_gtf:
    """
    Rename genes in combined TSEBRA GTF
    """
    input:
        rules.tsebra_combine.output
    output:
        gtf = f"{ANNOTATION_DIR}/braker/tsebra/braker_combined_renamed.gtf",
        tab = f"{ANNOTATION_DIR}/braker/tsebra/tsebra_rename_translationTab.txt"
    log: LOG_DIR + '/braker/rename_gtf.log'
    container: 'docker://teambraker/braker3:v.1.0.4'
    shell:
        """
        rename_gtf.py --gtf {input} \
                --prefix Psic \
                --translation_tab {output.tab} \
                --out {output.gtf} 2> {log} 
        """

###############################
###CLEAN & TRANSFORM GTF ####
###############################

rule remove_features:
    """
    Remove features except CDS. Required since BRAKER-generated GTFs are not compatible with most downstream software out-of-the box
    """
    input:
        rules.rename_tsebra_gtf.output
    output:
        f"{ANNOTATION_DIR}/cleaned/braker_combined_CDSonly_noOrgs.gtf"
    run:
        features = ['exon', 'intron', 'gene', 'transcript']
        with open(input[0], 'r') as fin:
            with open(output[0], 'w') as fout:
                lines = fin.readlines()
                for line in lines:
                    sline = line.split('\t')
                    if sline[2] in features:
                        pass
                    else:
                        fout.write(line)

rule gtf_to_gff:
    """
    Convert BRAKER-generated GTF to GFF3 using AGAT
    """
    input:
        rules.remove_features.output
    output:
        f"{ANNOTATION_DIR}/cleaned/Psic_structural.gff"
    container: 'docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0'
    log: LOG_DIR + '/gtf_to_gff/gtf_to_gff.log'
    shell:
        """
        agat_convert_sp_gxf2gxf.pl --gtf {input} \
            --output {output} &> {log}
        """

rule fix_seq:
    """
    Fix ID column where scaffold 0 was converted to SEQ
    """
    input:
        rules.gtf_to_gff.output
    output:
        f"{ANNOTATION_DIR}/cleaned/Psic_structural_seqFix.gff"
    shell:
        """
        sed 's/SEQ/0/g' {input} > {output}
        """

rule remove_internalStop_gene:
    """
    Remove single gene with internal stop codong causing InterproScan to break
    """
    input:
        rules.fix_seq.output
    output:
        f"{ANNOTATION_DIR}/cleaned/Psic_structural_seqFix_noInternalStop.gff"
    shell:
        """
        grep -v 'Psic_g38722' {input} > {output}
        """

rule gff_sort:
    """
    Sort GFF3 using genome tools
    """
    input:
        rules.remove_internalStop_gene.output
    output:
        f"{ANNOTATION_DIR}/cleaned/Psic_structural_sorted.gff"
    log: LOG_DIR + '/gff_sort/gff_sort.log'
    conda: '../envs/annotation.yaml'
    shell:
        """
        gt gff3 -sort -tidy -retainids {input} > {output} 2> {log}
        """

rule reformat_gff:
    """
    Fix transcript IDs for genes with alternative mRNA isoforms and remove transcript_id from gene features
    """
    input:
        rules.gff_sort.output
    output:
        f"{ANNOTATION_DIR}/cleaned/Psic_structural_sorted_reformated.gff"
    run:
        with open(input[0], 'r') as fin:
            with open(output[0], 'w') as fout:
                lines = fin.readlines()
                for line in lines:
                    if not line.startswith('#'):
                        sline = line.split('\t')
                        feature = sline[2]
                        if feature == 'gene':
                            #Remove transcript_id attribute from gene features
                            sline[8] = re.sub(r'(;transcript_id.*$)', '', sline[8])
                        elif feature == 'mRNA':
                            #Make sure transcript_id annotations for isoforms are correct (i.e., .t2, .t3, etc.)
                            #Use ID attribute since transcript_id attributes are incorrectly incremented and will be replaced
                            id_pattern = r"(?<=ID=)(.*)(?=;Parent)"
                            ID = re.search(id_pattern, sline[8]).group(1)
                            if ID.endswith('.t1'):
                                #First isoforms are fine
                                pass
                            else:
                                #Alternative isoforms need transcript_id replaced with ID
                                sline[8] = re.sub(r'(?<=;transcript_id=)(.*$)', ID, sline[8])
                        fout.write('\t'.join(sline))
                    else:
                        fout.write(line) 


rule get_proteins:
    """
    Get protein FASTA file using gffread
    """
    input:
        gff = rules.reformat_gff.output,
        ref = rules.repeat_masker.output.fasta
    output:
        f"{ANNOTATION_DIR}/Psic_proteins.fasta"
    log: LOG_DIR + '/get_proteins/get_proteins.log'
    container: 'docker://teambraker/braker3:v.1.0.4'
    shell:
        """
        gffread -E -y {output} -g {input.ref} {input.gff} 2> {log}
        """

###############################
###FUNCTIONAL ANNOTATION ####
###############################

rule run_interproscan:
    """
    Generate functional annotations using InterProScan
    """
    input:
        data = IPRSCAN_DATA,
        prot = rules.get_proteins.output
    output:
        gff =  f"{ANNOTATION_DIR}/interproscan/Psic_interproscan.gff3",
        xml =  f"{ANNOTATION_DIR}/interproscan/Psic_interproscan.xml"
    log: LOG_DIR + '/interproscan/run_interproscan.log'
    threads: 24
    params:
        out_base =  f"{ANNOTATION_DIR}/interproscan/Psic_interproscan"
    container: 'library://james-s-santangelo/interproscan/interproscan:5.61-93.0' 
    shell:
        """
        interproscan.sh -i {input.prot} \
            -b {params.out_base} \
            -f xml,gff3 \
            -goterms \
            --pathways \
            --seqtype p \
            --cpu {threads} \
            --verbose &> {log} 
        """ 

rule funannotate_setup:
    """
    Download Databases for Funannotate
    """
    output:
        directory(f"{ANNOTATION_DIR}/funannotate/fun_db")
    log: LOG_DIR + '/funannotate/funannotate_setup.log'
    container: 'docker://nextgenusfs/funannotate:v1.8.15'
    shell:
        """
        funannotate setup --database {output} -b tetrapoda --force 2> {log}
        """

rule dl_eggnog_db:
    """
    Download databases for EggNog-mapper
    """
    output:
        directory(f"{ANNOTATION_DIR}/eggnog/eggnog_db")
    conda: '../envs/annotation.yaml'
    shell:
        """
        mkdir -p {output}
        download_eggnog_data.py -y --data_dir {output} 
        """

rule run_eggnog_mapper:
    """
    Generate functional annotations using Eggnog-mapper
    """
    input:
        prot = rules.get_proteins.output,
        db = rules.dl_eggnog_db.output
    output:
        annot = f"{ANNOTATION_DIR}/eggnog/Psic.emapper.annotations"
    log: LOG_DIR + '/eggnog/aggnog_mapper.log'
    threads: 24
    conda: '../envs/annotation.yaml'
    params:
        out = f"{ANNOTATION_DIR}/eggnog/Psic"
    shell:
        """
        emapper.py -i {input.prot} \
            --data_dir {input.db} \
            --cpu {threads} \
            --itype proteins \
            -o {params.out} \
            --override &> {log}
        """

rule funannotate_annotate:
    """
    Use Funannotate to combine InterProScan and Eggnog annotations and generate additional functional annotations
    """
    input:
        enm = rules.run_eggnog_mapper.output.annot,
        gff = rules.reformat_gff.output, 
        ref = rules.repeat_masker.output.fasta,
        iprs = rules.run_interproscan.output.xml,
        db = rules.funannotate_setup.output 
    output:
        directory(f"{ANNOTATION_DIR}/funannotate/annotations"),
        fasta = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Podacris_siculus.scaffolds.fa", 
        prot = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Podacris_siculus.proteins.fa", 
        gff3 = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Podacris_siculus.gff3", 
        agp = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Podacris_siculus.agp", 
        gbk = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Podacris_siculus.gbk", 
        tbl = f"{ANNOTATION_DIR}/funannotate/annotations/annotate_results/Podacris_siculus.tbl"
    log: LOG_DIR + '/funannotate/funannotate_annotate.log'
    container: 'docker://nextgenusfs/funannotate:v1.8.15'
    threads: 24 
    params:
        outdir = f"{ANNOTATION_DIR}/funannotate/annotations"
    shell:
        """
        mkdir -p {params.outdir}
        funannotate annotate \
            --gff {input.gff} \
            --fasta {input.ref} \
            --species "Podacris siculus" \
            --out {params.outdir} \
            --iprscan {input.iprs} \
            --eggnog {input.enm} \
            --force \
            --database {input.db} \
            --busco_db tetrapoda \
            --cpus {threads} 2> {log}
        """


rule download_ec_numbers:
    """
    Download Enzyme Commission numbers from ExPASSY
    """
    output:
        f"{PROGRAM_RESOURCE_DIR}/EC_numbers/enzyme.dat"
    params:
        outdir = f"{PROGRAM_RESOURCE_DIR}/EC_numbers",
        url = 'https://ftp.expasy.org/databases/enzyme/enzyme.dat'
    shell:
        """
        wget {params.url} --no-check-certificate -P {params.outdir}
        """

rule fixEC_incrementCDS_addLocusTags:
    """
    Remap Hypothetical Proteins based on fully-resolved EC Numbers. Delete EC Number of not fully-resolved.
    Increment CDS IDs so they're unique.
    """
    input:
        gff = rules.funannotate_annotate.output.gff3,
        ec = rules.download_ec_numbers.output
    output:
        f"{ANNOTATION_DIR}/cleaned/Psic_functional_ECfix_wLocusTags.gff"
    conda: '../envs/annotation.yaml'
    params:
        locus_tag = 'Psic'
    script:
        "../scripts/python/fixEC_incrementCDS_addLocusTags.py"

rule gff_sort_functional:
    """
    Sort GFF3 with functional annotations using GFF3_sort Perl script. Can't use genometools here since sorting not compatible with table2asn
    """
    input:
        rules.fixEC_incrementCDS_addLocusTags.output
    output:
        f"{ANNOTATION_DIR}/Psic_functional_final_sorted.gff3"
    log: LOG_DIR + '/gff_sort/gff_sort_functional.log'
    conda: '../envs/annotation.yaml'
    shell:
        """
        gff3_sort -g {input} -r -og {output} 2> {log}
        """

rule get_proteins_finalGFF:
    """
    Get final protein set using GFFREAD
    """
    input:
        gff = rules.gff_sort_functional.output,
        ref = rules.repeat_masker.output.fasta
    output:
        f"{ANNOTATION_DIR}/Psic_proteins_final.fasta"
    log: LOG_DIR + '/get_proteins/get_proteins_final.log'
    container: 'docker://teambraker/braker3:v.1.0.4'
    shell:
        """
        gffread -E -y {output} -g {input.ref} {input.gff} 2> {log}
        """

rule annotation_done:
    input:
        rules.get_proteins_finalGFF.output
    output:
        f"{ANNOTATION_DIR}/annotation.done"
    shell:
        """
        touch {output}
        """
