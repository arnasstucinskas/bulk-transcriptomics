SAMPLES = [
    "KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-3_S8",
    "KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-2_S7",
    "KAPA_mRNA_HyperPrep_-HBR-KAPA-100_ng_total_RNA-3_S6",
    "KAPA_mRNA_HyperPrep_-HBR-KAPA-100_ng_total_RNA-2_S5",
    "Collibri_standard_protocol-UHRR-Collibri-100_ng-2_S3",
    "Collibri_standard_protocol-UHRR-Collibri-100_ng-3_S4",
    "Collibri_standard_protocol-HBR-Collibri-100_ng-3_S2",
    "Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1"
]

rule all:
    input:
        expand("outputs/fastqc_results/{sample}_L001_R1_001_fastqc.html", sample=SAMPLES),
        expand("outputs/fastqc_results/{sample}_L001_R2_001_fastqc.html", sample=SAMPLES),
        expand("outputs/trimmed_data/{sample}_L001_R1_001_trimmed.fastq", sample=SAMPLES),
        expand("outputs/trimmed_data/{sample}_L001_R2_001_trimmed.fastq", sample=SAMPLES),
        expand("outputs/trimmed_fastqc_results/{sample}_L001_R1_001_trimmed_fastqc.html", sample=SAMPLES),
        expand("outputs/trimmed_fastqc_results/{sample}_L001_R2_001_trimmed_fastqc.html", sample=SAMPLES),
        "outputs/multiqc_report.html",
        "outputs/trimmed_multiqc_report.html",
        expand("outputs/aligned_reads/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("outputs/aligned_reads/{sample}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand("outputs/counts/{sample}_counts_strand1.txt", sample=SAMPLES),
        expand("outputs/counts/{sample}_counts_strand2.txt", sample=SAMPLES),
        expand("outputs/counts/{sample}_best_strand_number.txt", sample=SAMPLES),
        expand("outputs/counts/{sample}_best_strand.txt", sample=SAMPLES),
        "outputs/deseq2/count_matrix.csv",
        "outputs/deseq2/column_data.csv",
        "outputs/deseq2/deseq2_results.csv",
        "outputs/deseq2/volcano_plot.png",
        "outputs/deseq2/pca_plot_de_genes.png",
        "outputs/deseq2/enrichment_results.csv"
        

#################################################################################################
# 1. Evaluate  "raw" reads quality via FastQC.
#################################################################################################

# FastQC analysis rule
rule fastqc:
    input:
        r1="data/{sample}_L001_R1_001.fastq.gz",
        r2="data/{sample}_L001_R2_001.fastq.gz"
    output:
        r1_html="outputs/fastqc_results/{sample}_L001_R1_001_fastqc.html",
        r2_html="outputs/fastqc_results/{sample}_L001_R2_001_fastqc.html"
    shell:
        "fastqc {input.r1} {input.r2} --outdir ./outputs/fastqc_results"

#################################################################################################
# 2. Create a raport via MultiQC.
#################################################################################################

# Rule for generating a MultiQC report
rule multiqc:
    input:
        expand("outputs/fastqc_results/{sample}_L001_R1_001_fastqc.html", sample=SAMPLES),
        expand("outputs/fastqc_results/{sample}_L001_R2_001_fastqc.html", sample=SAMPLES)
    output:
        "outputs/multiqc_report.html"
    shell:
        "multiqc outputs/fastqc_results/ --outdir ./outputs"

#################################################################################################
# 3. Trim a barcodes with bbduk.
#################################################################################################

# Rule for trimming barcodes with bbduk
rule trim_barcodes:
    input:
        r1="data/{sample}_L001_R1_001.fastq.gz",
        r2="data/{sample}_L001_R2_001.fastq.gz"
    output:
        trimmed_r1="outputs/trimmed_data/{sample}_L001_R1_001_trimmed.fastq",
        trimmed_r2="outputs/trimmed_data/{sample}_L001_R2_001_trimmed.fastq"
    params:
        adapters="adapters.fa",
        ktrim="r",
        k=23,
        mink=11,
        hdist=1,
        tpe=True,
        tbo=True,
        qtrim="r",
        trimq=10
    shell:
        "bbduk.sh in1={input.r1} in2={input.r2} "
        "out1={output.trimmed_r1} out2={output.trimmed_r2} "
        "ref={params.adapters} ktrim={params.ktrim} k={params.k} "
        "mink={params.mink} hdist={params.hdist} tpe={params.tpe} tbo={params.tbo} "
        "qtrim={params.qtrim} trimq={params.trimq}"

#################################################################################################
# 4. Run FastQC analysis for the trimmed reads.
#################################################################################################

# Rule for FastQC analysis on trimmed reads
rule fastqc_trimmed:
    input:
        r1="outputs/trimmed_data/{sample}_L001_R1_001_trimmed.fastq",
        r2="outputs/trimmed_data/{sample}_L001_R2_001_trimmed.fastq"
    output:
        r1_html="outputs/trimmed_fastqc_results/{sample}_L001_R1_001_trimmed_fastqc.html",
        r2_html="outputs/trimmed_fastqc_results/{sample}_L001_R2_001_trimmed_fastqc.html"
    shell:
        "fastqc {input.r1} {input.r2} --outdir outputs/trimmed_fastqc_results"

#################################################################################################
# 5. Create a separate raport on the trimmed reads raport for MultiQC.
# Compare FastQC/MultiQC output for trimmed and not trimmed reads.
# What indicates that adapters were successfully removed.
#################################################################################################

rule multiqc_trimmed:
    input:
        expand("outputs/trimmed_fastqc_results/{sample}_L001_R1_001_trimmed_fastqc.html", sample=SAMPLES),
        expand("outputs/trimmed_fastqc_results/{sample}_L001_R2_001_trimmed_fastqc.html", sample=SAMPLES)
    output:
        "outputs/trimmed_multiqc_report.html"
    shell:
        "multiqc outputs/trimmed_fastqc_results/ --outdir ./outputs -n trimmed_multiqc_report.html"

#################################################################################################
# 6. Add STAR aligner step to the workflow (this would include two steps and genome indexing).
#################################################################################################

# a) Rule for indexing the genome with STAR (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
rule star_index:
    input:
        fa="data/play_data_ref_annot/chr19_20Mb.fa",
        gtf="data/play_data_ref_annot/chr19_20Mb.gtf"
    output:
        directory("outputs/star_index")
    shell:
        "STAR "
        "--runMode genomeGenerate "
        "--genomeSAindexNbases 11 "
        "--runThreadN 4 "
        "--genomeDir outputs/star_index "
        "--genomeFastaFiles data/play_data_ref_annot/chr19_20Mb.fa "
        "--sjdbGTFfile data/play_data_ref_annot/chr19_20Mb.gtf"

# b) Rule for aligning trimmed reads with STAR (use --outSAMtype BAM SortedByCoordinate to get a sorted alignment file)
rule star_align:
    input:
        index="outputs/star_index",
        r1="outputs/trimmed_data/{sample}_L001_R1_001_trimmed.fastq",
        r2="outputs/trimmed_data/{sample}_L001_R2_001_trimmed.fastq"
    output:
        bam="outputs/aligned_reads/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        prefix="outputs/aligned_reads/{sample}_"
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {input.index} "
        "--readFilesIn {input.r1} {input.r2} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM SortedByCoordinate"

# #################################################################################################
# # 7. Add a step that indexes the output bam with samtools index
# # Add counting step with featureCounts http://bioinf.wehi.edu.au/featureCounts/
# #################################################################################################

# a) Rule for indexing BAM files with samtools
rule samtools_index:
    input:
        "outputs/aligned_reads/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "outputs/aligned_reads/{sample}_Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

# B) Rule for read counting with featureCounts
rule feature_counts:
    input:
        bam="outputs/aligned_reads/{sample}_Aligned.sortedByCoord.out.bam",
        gtf="data/play_data_ref_annot/chr19_20Mb.gtf"
    output:
        "outputs/counts/{sample}_counts_strand{strand}.txt"
    params:
        strand_specificity=lambda wildcards: wildcards.strand
    wildcard_constraints:
        strand="1|2"
    shell:
        "featureCounts -p -t exon -g gene_id "
        "-a {input.gtf} -o {output} {input.bam} -s {params.strand_specificity}"

# Select best strand number
rule compare_feature_counts:
    input:
        counts1_summary="outputs/counts/{sample}_counts_strand1.txt",
        counts2_summary="outputs/counts/{sample}_counts_strand2.txt"
    output:
        best_strand_number_file="outputs/counts/{sample}_best_strand_number.txt",
        best_strand_file="outputs/counts/{sample}_best_strand.txt"
    script:
        "scripts/compare_counts.py"

rule prepare_deseq2_input:
    input:
        counts=expand("outputs/counts/{sample}_best_strand.txt", sample=SAMPLES),
    output:
        count_matrix="outputs/deseq2/count_matrix.csv",
        col_data="outputs/deseq2/column_data.csv"
    script:
        "scripts/prepare_deseq2_input.py"

# #################################################################################################
# # 8. For each sample preparation method perform a DE analysis comparing UHRRvsHBR. 
# # UHRR is a sample originating from several cancerous cell lines,HBR represents “normal” sample
# #################################################################################################

# 8.1. DESeq2 differential expression analysis
rule deseq2_analysis:
    input:
        count_matrix="outputs/deseq2/count_matrix.csv",
        col_data="outputs/deseq2/column_data.csv"
    output:
        deseq2_results="outputs/deseq2/deseq2_results.csv"
    script:
        "scripts/run_deseq2_analysis.R"

# 8.2. Generating a volcano plot from the DESeq2 results
rule volcano_plot:
    input:
        deseq2_results="outputs/deseq2/deseq2_results.csv"
    output:
        volcano_plot="outputs/deseq2/volcano_plot.png"
    script:
        "scripts/generate_volcano_plot.R"

# #################################################################################################
# # 9. Create a PCA plot exploring how different sample preparation plots cluster based on 
# # deferentially expressed genes (use genes that are DE for all sample preparation methods)
# #################################################################################################
rule pca_plot_de_genes:
    input:
        deseq2_results="outputs/deseq2/deseq2_results.csv"
    output:
        pca_plot="outputs/deseq2/pca_plot_de_genes.png"
    script:
        "scripts/generate_pca_plot.R"

# #################################################################################################
# # 10. Pathway enrichment analysis using fgsea
# #################################################################################################
rule fgsea_analysis:
    input:
        deseq2_results="outputs/deseq2/deseq2_results.csv"
    output:
        enrichment_results="outputs/deseq2/enrichment_results.csv",
        fgsea_analysis_plot="outputs/deseq2/fgsea_analysis_plot.png"
    script:
        "scripts/run_fgsea_analysis.R"

