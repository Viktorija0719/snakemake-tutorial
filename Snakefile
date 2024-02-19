configfile: "config.yaml"

rule all:
    input:
        "data/complete",
        "plots/quals.svg"


rule download_data:
    output:
        "data/snakemake-tutorial-data.tar.gz"
    shell:
        "wget -O {output} https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz"

rule extract_data:
    input:
        "data/snakemake-tutorial-data.tar.gz"
    output:
        touch("data/complete") 
    shell:
        """
        tar -xzvf {input} -C data/
        mv data/snakemake-tutorial-data-*/data/* data/
        rm -r data/snakemake-tutorial-data-*
        touch {output}
        """

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule bwa_map:
    input:
        fa="data/genome.fa",
        fq=get_bwa_map_input_fastqs,
        complete="data/complete"  # Ensure this rule waits for data preparation
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input.fa} {input.fq} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    params:
        rate=config["prior_mutation_rate"]
    log:
        "logs/bcftools_call/all.log"
    shell:
        "(bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv -P {params.rate} - > {output}) 2> {log}"

rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
