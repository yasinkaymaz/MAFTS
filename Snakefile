
dockerTag = "latest" #FIXME tagged versions

rule all:
    input:
        expand("FastQC/{sample}_R{num}_fastqc.zip", sample=config["samples"], num=['1', '2']),
        #"utdir/fastq_multiqc.html",
        expand("gatk/{sample}_SnpEff.vcf.gz", sample=config["samples"]),
        expand("gatk/{sample}_snpEff_summary.html", sample=config["samples"])

rule FastQC:
    input:
        fq=lambda wildcards: f"{config['samples'][wildcards.sample]}_R{wildcards.num}.fastq.gz"
    output:
        html="FastQC/{sample}_R{num}_fastqc.html",
        zip="FastQC/{sample}_R{num}_fastqc.zip"
    shell:
        "fastqc "
        "{input.fq} "
        "--outdir FastQC"

rule MultiFastQC:
    input:
        expand("FastQC/{sample}_R{num}_fastqc.html", sample=config["samples"], num=['1', '2'])
    output:
        "utdir/fastq_multiqc.html"
    shell:
        "multiqc FastQC/ -outdir utdir --filename fastq_multiqc.html"

rule AdpTrimFastp:
    input:
        sample=["fastq_files/{sample}_R1.fastq.gz", "fastq_files/{sample}_R2.fastq.gz"]
    output:
        trimmed=["trimmed/{sample}_R1.fastq", "trimmed/{sample}_R2.fastq"],
        unpaired1=temp("trimmed/{sample}.u1.fastq"),
        unpaired2=temp("trimmed/{sample}.u2.fastq"),
        merged=temp("trimmed/{sample}.merged.fastq"),
        failed=temp("trimmed/{sample}.failed.fastq"),
        html="report/{sample}.html",
        json="report/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    threads: 2
    wrapper:
        "v0.86.0/bio/fastp"

rule Bowtie2:
    input:
        sample=["trimmed/{sample}_R1.fastq", "trimmed/{sample}_R2.fastq"]
    output:
        temp("mapped/{sample}.bam")
    log:
        "logs/bowtie2/{sample}.log"
    params:
        index=config["GenomeIndex"],  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: 6  # Use at least two threads
    wrapper:
        "v0.86.0/bio/bowtie2/align"

rule Samtools_sort:
    input:
        "mapped/{sample}.bam"
    output:
        temp("mapped/{sample}.sorted.bam")
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads:
        4
    wrapper:
        "v0.86.0/bio/samtools/sort"

rule chrFilter:
    input:
        ibam="mapped/{sample}.sorted.bam"
    output:
        file=temp("mapped/{sample}.chr1_22.bam")
    shell:
        "samtools index {input.ibam} && "
        "samtools view -o {output.file} {input.ibam} `seq 1 22 | sed 's/^/chr/'` "

rule Samtools_coverage:
    input:
        bam="mapped/{sample}.chr1_22.bam",
        bed=config["TargetRegions"]
    output:
        cov="samtools/cov/{sample}.coverage.txt"
    params:
        extra="" # optional additional parameters as string
    shell:
        "bedtools coverage -a {input.bed} -b {input.bam} > {output.cov} "

rule BamAddReadGroups:
    input:
        bam="mapped/{sample}.chr1_22.bam"
    output:
        outbam="mapped/{sample}.RG.bam"
#    container:
#        "docker://broadinstitute/gatk:{}".format(dockerTag)
    threads: 4
    log: 'logs/gatk/mutect2-tumor-only/{sample}.log'
    params:
        readgroup='RG-{sample}',
        gatkDir=config["GATKdir"]
    shell:
        "{params.gatkDir}/gatk AddOrReplaceReadGroups "
        "-I {input.bam} "
        "-O {output.outbam} "
        "-LB {params.readgroup} "
        "-PU {params.readgroup} "
        "-SM {params.readgroup} "
        "-PL 'ILLUMINA' "

rule Mutect2_tumor_only:
    input:
        bam="mapped/{sample}.RG.bam",
        genome=config["GenomeFa"]
    output:
        "gatk/{sample}.vcf.gz"
    params:
        gatkDir=config["GATKdir"]
    shell:
        "samtools index {input.bam} && "
        "{params.gatkDir}/gatk Mutect2 "
        "-I {input.bam} "
        "-R {input.genome} "
        "-O {output} "

rule FilterCalls:
    input:
        genome=config["GenomeFa"],
        vcf="gatk/{sample}.vcf.gz"
    output:
        temp("gatk/{sample}_filtered.vcf.gz")
    params:
        gatkDir=config["GATKdir"],
        tabDir=config["tabixDir"]
    shell:
        "{params.gatkDir}/gatk VariantFiltration "
        "-V {input.vcf} "
        "-R {input.genome} "
        "--filter-expression  'DP <= 20' --filter-name 'DepthofQuality' "
        "-O {output}.tmp && "
        "bcftools view -f PASS {output}.tmp |{params.tabDir}bgzip > {output} && "
        "bcftools index --tbi {output} "

rule FilterMutectCalls:
    input:
        genome=config["GenomeFa"],
        vcf="gatk/{sample}_filtered.vcf.gz"
    output:
        temp("gatk/{sample}_MutFiltered.vcf.gz")
    params:
        statfile="gatk/{sample}.vcf.gz.stats",
        gatkDir=config["GATKdir"],
        tabDir=config["tabixDir"]
    shell:
        "{params.gatkDir}/gatk FilterMutectCalls "
        "-V {input.vcf} "
        "-R {input.genome} "
        "--stats {params.statfile} "
        "--min-allele-fraction 0.1 "
        "--unique-alt-read-count 5 "
        "-O {output}.tmp && "
        "bcftools view -f PASS {output}.tmp |{params.tabDir}bgzip > {output} && "
        "bcftools index --tbi {output} "

rule AnnotateSNPs:
    input:
        vcf="gatk/{sample}_MutFiltered.vcf.gz",
        bam="mapped/{sample}.RG.bam"
    output:
        temp("gatk/{sample}_SNPFiltered.vcf.gz")
    params:
        dbsnpFile=config["dbSNPfile"],
        gatkDir=config["GATKdir"]
    shell:
        "{params.gatkDir}/gatk VariantAnnotator "
        "-V {input.vcf} "
        "--dbsnp {params.dbsnpFile} "
        "-O {output} "

rule SnpEff:
    input:
        vcf="gatk/{sample}_SNPFiltered.vcf.gz"
    output:
        vcEff="gatk/{sample}_SnpEff.vcf.gz",
        genesTxt="gatk/{sample}_snpEff_genes.txt",
        summary="gatk/{sample}_snpEff_summary.html"
    params:
        tabDir=config["tabixDir"]
    shell:
        "snpEff GRCh38.99 {input.vcf}|{params.tabDir}bgzip > {output.vcEff} && "
        "mv snpEff_genes.txt {output.genesTxt} && "
        "mv snpEff_summary.html {output.summary} "
