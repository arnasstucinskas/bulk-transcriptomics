SAMPLES = [
    # "KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-3_S8",
    # "KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-2_S7",
    # "KAPA_mRNA_HyperPrep_-HBR-KAPA-100_ng_total_RNA-3_S6",
    # "KAPA_mRNA_HyperPrep_-HBR-KAPA-100_ng_total_RNA-2_S5",
    # "Collibri_standard_protocol-UHRR-Collibri-100_ng-2_S3",
    # "Collibri_standard_protocol-UHRR-Collibri-100_ng-3_S4",
    # "Collibri_standard_protocol-HBR-Collibri-100_ng-3_S2",
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
        expand("outputs/aligned_reads/{sample}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES), # step 7
        expand("outputs/counts/{sample}_counts.txt", sample=SAMPLES)

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
        fa="data/play_data_ref_annot/chr19_20Mb.fa",  # Update this path to your genome fasta file
        gtf="data/play_data_ref_annot/chr19_20Mb.gtf"  # Update this path to your genome GTF file
    output:
        directory("outputs/star_index")
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir outputs/star_index "
        "--genomeFastaFiles {input.fa} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99"  # Set this to the length of your reads minus 1

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
        "--readFilesCommand zcat "
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
        "outputs/counts/{sample}_counts.txt"
    params:
        strand_specificity="2"  # 1 or 2 based on library
    shell:
        "featureCounts -p -t exon -g gene_id "
        "-a {input.gtf} -o {output} {input.bam} -s {params.strand_specificity}"

