#single-end

import os
import pandas as pd
from os import path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

configfile: "envs/config.yaml"

STUDY = config["study6"]["name"]
LOG_DIR = "logs/" + STUDY + "/"
WORK_DIR = "data/" + STUDY + "/"
RESULT_DIR = "results/" + STUDY + "/"
BENCHMARK_DIR="benchmarks/" + STUDY + "/"
READ_LEN = config["study6"]["length"]
df = pd.read_csv(config["study6"]["route"], sep="\t")
samples = df["Run"].tolist()

LEVEL = ["G", "S"]
TYPE = ["bacteria", "virus"]

def get_level(level):
	return pd.DataFrame({"values": [0,1], "level": level})


rule all:
	input:
#		expand(RESULT_DIR + "fastq/{sample}.sra.fastq",sample=samples),
#		expand(RESULT_DIR + "trimmed/{sample}.{a}",a=["sra_trimmed.fq.gz","sra.fastq_trimming_report.txt"],sample=samples),
#		expand(RESULT_DIR + "fastqc/{sample}.sra_fastqc.{ext}",ext=["html","zip"], sample=samples),
#		expand(RESULT_DIR + "trimmed_qc/{sample}.sra_trimmed_fastqc.{ext}",ext=["html","zip"], sample=samples),
		RESULT_DIR + "multiqc/multiqc.html",
		RESULT_DIR + "multiqc/trim_multiqc.html",
#		expand(RESULT_DIR + "star_align/{sample}_", sample=samples),
		expand(RESULT_DIR + "kraken2_classify/{sample}_{type}.{ext}", ext=["kreport","output"], sample=samples, type=TYPE),
#		expand(RESULT_DIR + "markduplicate/{sample}_marked_duplicates.bam",sample=samples),
#		expand(RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam",sample=samples),
#		expand(RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam.bai",sample=samples),
#		expand(RESULT_DIR+"htseq_count/{sample}.count",sample=samples),
		RESULT_DIR+"rRNA_percentage/rRNA_rate.csv",
		expand(RESULT_DIR + "bracken_estimation/{sample}_{type}_{level}.{ext}", ext=["bracken","breport","mpa.txt"], sample=samples, level=LEVEL, type=TYPE)
		expand(RESULT_DIR + "combine_result/{type}_{level}.{ext}", ext=["bracken.txt","mpa.txt"], level=LEVEL, type=TYPE),
		expand(RESULT_DIR + "htseq_count_exon/{sample}.count",sample=samples),
		expand(RESULT_DIR + "read_distribute/{sample}_read_distribute.log",sample=samples)
		


rule Fasterq_dump:
	input:
		WORK_DIR + "{sample}/{sample}.sra"
	output:
		temp(RESULT_DIR + "fastq/{sample}.sra.fastq")
	params:
		outdir =RESULT_DIR + "fastq"
	threads: 2
	benchmark:
		BENCHMARK_DIR + "fastq/{sample}.fasterq_dump.benchmark.txt"
	log:
		LOG_DIR + "fastq/{sample}.sra2fq.log"
	conda:
		"envs/fastq.yaml"
	priority: 10
	message:
		"Executing fasterq-dump with {threads} threads on the following files {input}."
	shell:
		"""
		(fasterq-dump --split-3 -e {threads} {input} -o {output}) &> {log}
		"""
rule TrimGalore:
	input:
		RESULT_DIR + "fastq/{sample}.sra.fastq"
	output:
		temp(RESULT_DIR + "trimmed/{sample}.sra_trimmed.fq.gz"),
		RESULT_DIR + "trimmed/{sample}.sra.fastq_trimming_report.txt"
	params:
		outdir = RESULT_DIR + "trimmed/"
	benchmark:
                BENCHMARK_DIR + "trimmed/{sample}.trimmed.benchmark.txt"
	conda:
		"envs/trimgalore.yaml"
	log:
		LOG_DIR + "trim_galore/{sample}.trimmed.log"
	message:
		"Executing TrimGalore on the following files {input}."
	priority: 9
	shell:
		"(trim_galore --phred33 --gzip --output_dir {params.outdir} {input}) &> {log}"
rule FastQC:
	input:
		fq = RESULT_DIR + "fastq/{sample}.sra.fastq",
		trim_fq = RESULT_DIR + "trimmed/{sample}.sra_trimmed.fq.gz"
	output:
		RESULT_DIR + "fastqc/{sample}.sra_fastqc.html",
		RESULT_DIR + "fastqc/{sample}.sra_fastqc.zip",
		RESULT_DIR + "trimmed_qc/{sample}.sra_trimmed_fastqc.html",
		RESULT_DIR + "trimmed_qc/{sample}.sra_trimmed_fastqc.zip"
	benchmark:
		BENCHMARK_DIR+"fastqc/{sample}.fastqc.benchmark.txt",
	conda:
		"envs/fastqc.yaml"
	params:
		outdir1 = RESULT_DIR + "fastqc",
		outdir2 = RESULT_DIR + "trimmed_qc"
	log:
		log1 = LOG_DIR + "fastqc/{sample}.fastqc.log",
		log2 = LOG_DIR + "trimmed_qc/{sample}.trim_qc.log"
	message:
		"Executing fastqc on the following files {input}."
	priority: 8
	shell:
		"""
		fastqc -o {params.outdir1} --noextract {input.fq} &> {log.log1}
		fastqc -o {params.outdir2} --noextract {input.trim_fq} &> {log.log2}
		"""
rule MultiQC:
	input:
		expand(RESULT_DIR + "fastqc/{sample}.sra_fastqc.{ext}",ext=["html","zip"], sample=samples),
		expand(RESULT_DIR + "trimmed_qc/{sample}.sra_trimmed_fastqc.{ext}",ext=["html","zip"], sample=samples)
	output:
		RESULT_DIR + "multiqc/multiqc.html",
		RESULT_DIR + "multiqc/trim_multiqc.html"
	params:
		outdir = RESULT_DIR + "multiqc",
		indir1 = RESULT_DIR + "fastqc",
		outname1 = "multiqc.html",
		indir2 = RESULT_DIR + "trimmed_qc",
		outname2 = "trim_multiqc.html"
	conda:
		"envs/multiqc.yaml"
	log:
		log1 = LOG_DIR + "multiqc/multiqc.log",
		log2 = LOG_DIR + "multiqc/trim_multiqc.log"
	message:
		"Executing MultiQC"
	priority: 7
	shell:
		"""
		multiqc --force -o {params.outdir} -n {params.outname1} {params.indir1} &> {log.log1}
		multiqc --force -o {params.outdir} -n {params.outname2} {params.indir2} &> {log.log2}
		"""
rule STAR_Align:
	input:
		RESULT_DIR + "trimmed/{sample}.sra_trimmed.fq.gz"
	output:
		RESULT_DIR + "star_align/{sample}_Aligned.sortedByCoord.out.bam",
		RESULT_DIR + "star_align/{sample}_Aligned.toTranscriptome.out.bam",
		RESULT_DIR + "star_align/{sample}_Unmapped.out.mate1",
		RESULT_DIR + "star_align/{sample}_"
	params:
		index = config["STAR_index"],
		sample = RESULT_DIR + "star_align/{sample}_"
	benchmark:
		BENCHMARK_DIR + "star_align/{sample}.star_align.benchmark.txt"
	threads: 12
	resources:
		mem_mb = 30000
	conda:
		"envs/star_align.yaml"
	log:
		LOG_DIR + "star/{sample}_star_align.log"
	message:
		"Executing STAR Align on the following files {input}."
	priority: 10
	shell:
		"""
		(STAR --runThreadN {threads} --runMode alignReads --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts \
		--twopassMode Basic --outSAMtype BAM SortedByCoordinate --outReadsUnmapped  Fastx --outBAMsortingThreadN 5 \
		--genomeDir {params.index}  --readFilesIn {input[0]} --outFileNamePrefix {params.sample}) &> {log}
		"""

rule Markduplicate:
	input:
		RESULT_DIR + "star_align/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		RESULT_DIR + "markduplicate/{sample}_marked_duplicates.bam",
		RESULT_DIR + "markduplicate/{sample}_marked_dup_metrics.txt"
	benchmark:
		BENCHMARK_DIR + "markduplicate/{sample}.markduplicate.benchmark.txt"
	log:
		LOG_DIR + "markduplicate/{sample}_md.log"
	conda:
		"envs/markduplicate.yaml"
	threads:2
	message:
		"Executing markdupliacte on the following files {input}."
	priority: 11
	shell:
		"(picard MarkDuplicates -I {input} -O {output[0]} -M {output[1]}) &> {log}"



rule Reads_distribution:
        input:
               RESULT_DIR + "markduplicate/{sample}_marked_duplicates.bam"
        output:
                RESULT_DIR + "read_distribute/{sample}_read_distribute.log"
        params:
                bed = config["bed"]
        benchmark:
                BENCHMARK_DIR + "read_distribute/{sample}.read_distribute.benchmark.txt"
        threads: 12
        conda:
                "envs/RseQC.yaml"
        message:
                "Executing read distribution on the following files {input}."
        priority: 13
        shell:
                "read_distribution.py -i {input} -r {params.bed} &> {output}"

rule BAM_filter:
	input:
		 RESULT_DIR + "markduplicate/{sample}_marked_duplicates.bam"
	output:
		RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam"
	benchmark:
                BENCHMARK_DIR + "markduplicate/{sample}.bamfilter.benchmark.txt"
	log:
		LOG_DIR+"BAM_filter/{sample}_BAM_filter.log"
	conda:
		"envs/bam_filter.yaml"
	threads:2
	message:
		"Executing BAM filter on the following files {input}."
	priority: 12
	shell:
		"(samtools view -f 0 -F 1024 {input} -o {output}) &> {log}"

rule Samtools_index:
	input:
		RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam"
	output:
		RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam.bai"
	benchmark:
		BENCHMARK_DIR + "BAM_filter/{sample}.samtools_index.benchmark.txt"
	log:
		LOG_DIR + "BAM_filter/{sample}_BAM_filter.log"
	conda:
		"envs/bam_filter.yaml"
	threads:2
	message:
		"Executing samtools index on the following files {input}."
	priority: 13
	shell:
		"(samtools index {input} > {output}) &> {log}"

rule Reads_count:
	input:
		bam = RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam",
		bai = RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam.bai",
		gtf=config["ref_gtf"]
	output:
		RESULT_DIR + "htseq_count/{sample}.count"
	benchmark:
               BENCHMARK_DIR + "htseq_count/{sample}.count.benchmark.txt"
	log:
		LOG_DIR + "htseq_count/{sample}.count.log"
	conda:
		"envs/htseq_count.yaml"
	message:
		"Executing reads count on the following files {input}."
	priority: 14
	shell:
		"(htseq-count  -f bam -r pos -s no -a 10 -t exon -i gene_id -m union --nonunique=none {input.bam} {input.gtf} > {output}) &>{log}"
#		"(htseq-count  -f bam -r pos -s yes -a 10 -t exon -i gene_id -i exon_number -m union --nonunique=none --additional-attr=gene_name --additional-attr=exon_number {input.bam} {input.gtf} > {output}) &>{log}"
#		"(htseq-count  -f bam -r pos -s reverse -a 10 -t exon -i gene_id -i exon_number -m union --nonunique=none --additional-attr=gene_name --additional-attr=exon_number {input.bam} {input.gtf} > {output}) &>{log}"


rule Exon_number:
	input:
		bam = RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam",
		bai = RESULT_DIR + "BAM_filter/{sample}_BAM_filter.bam.bai",
		gtf=config["ref_gtf"]
	output:
		RESULT_DIR + "htseq_count_exon/{sample}.count"
	benchmark:
		BENCHMARK_DIR + "htseq_count_exon/{sample}.count.benchmark.txt"
	log:
		LOG_DIR + "htseq_count_exon/{sample}.count.log"
	conda:
		"envs/htseq_count.yaml"
	message:
		"Executing reads count on the following files {input}."
	priority: 14
	shell:
		"(htseq-count  -f bam -r pos -s no -a 10 -t exon -i gene_id -i exon_number -m union --nonunique=none --additional-attr=gene_name --additional-attr=exon_number {input.bam} {input.gtf} > {output} ) &>{log}"
#		"(htseq-count  -f bam -r pos -s yes -a 10 -t exon -i gene_id -i exon_number -m union --nonunique=none --additional-attr=gene_name --additional-attr=exon_number {input.bam} {input.gtf} > {output}) &>{log}"
#		"(htseq-count  -f bam -r pos -s reverse -a 10 -t exon -i gene_id -i exon_number -m union --nonunique=none --additional-attr=gene_name --additional-attr=exon_number {input.bam} {input.gtf} > {output}) &>{log}"


rule RNA_degrade:
	input:
		expand(RESULT_DIR + "htseq_count_exon/{sample}.count",sample=samples)
	output:
		RESULT_DIR + "RNA_degradation/degrate_genes.csv"
	script:
		"scripts/RNA_degradation.R"

rule rRNA_percentage:
	input:
		expand(RESULT_DIR + "htseq_count/{sample}.count",sample=samples)
	output:
		RESULT_DIR + "rRNA_percentage/rRNA_rate.csv"
	script:
		"scripts/rRNA_percentage.R"

rule Kraken2_classify:
	input:
		fq = RESULT_DIR + "star_align/{sample}_Unmapped.out.mate1"
	output:
		krp = RESULT_DIR + "kraken2_classify/{sample}_{type}.kreport",
		kop = RESULT_DIR + "kraken2_classify/{sample}_{type}.output"
	benchmark:
                BENCHMARK_DIR + "kraken2_classify/{sample}_{type}.kraken2_classify.benchmark.txt"
	params:
		db = lambda w: config["kraken_database"][w.type]
	threads: 12
	resources:
		mem_mb = 45000
	conda:
		"envs/kraken2.yaml"
	log:
		log1 = LOG_DIR + "kraken2_classify/{sample}_bacteria_kraken2.log",
		log2 = LOG_DIR + "kraken2_classify/{sample}_viral_kraken2.log"
	message:
		"Executing Kraken2 on the following files {input}."
	priority: 3
	shell:
		"""
		kraken2 --db {params.db} --threads {threads} --report {output.krp} --output {output.kop} {input.fq}  &> {log}
		"""


rule Bracken:
	input:
		RESULT_DIR + "kraken2_classify/{sample}_{type}.kreport"
	output:
		brk = RESULT_DIR + "bracken_estimation/{sample}_{type}_{level}.bracken",
		brp = RESULT_DIR + "bracken_estimation/{sample}_{type}_{level}.breport",
		mpa = RESULT_DIR + "bracken_estimation/{sample}_{type}_{level}.mpa.txt"
	benchmark:
		BENCHMARK_DIR + "bracken_estimation/{sample}_{type}_{level}.bracken_estimation.benchmark.txt"
	params:
		db = lambda w: config["kraken_database"][w.type],
		level = lambda wc: wc.get("level"),
		readlen = READ_LEN
	threads:1
	conda:
		"envs/bracken.yaml"
	log:
		log1 = LOG_DIR + "bracken_estimation/{sample}_{type}_{level}_bracken.log",
		log2 = LOG_DIR + "bracken_estimation/{sample}_{type}_{level}_brp2mpa.log"
	message:
		"Executing Bracken on the following files {input}."
	priority: 1
	shell:
		"""
		bracken -d {params.db} -i {input} -o {output.brk} -w {output.brp} -l {params.level} -r {params.readlen} -t {threads} &> {log.log1}
		kreport2mpa.py --display-header -r {output.brp} -o {output.mpa} &> {log.log2}
		"""
rule Combine:
	input:
		brp = expand(RESULT_DIR + "bracken_estimation/{sample}_{type}_{level}.bracken",sample=samples, level=LEVEL, type=TYPE),
		mpa = expand(RESULT_DIR + "bracken_estimation/{sample}_{type}_{level}.mpa.txt",sample=samples, level=LEVEL, type=TYPE)
	output:
		op1 = RESULT_DIR + "combine_result/{type}_{level}.bracken.txt",
		op2 = RESULT_DIR + "combine_result/{type}_{level}.mpa.txt"
	benchmark:
		BENCHMARK_DIR + "combine_result/{type}_{level}.combine_result.benchmark.txt"
	params:
		input1 = RESULT_DIR + "bracken_estimation/*_{type}_{level}.bracken",
		input2 = RESULT_DIR + "bracken_estimation/*_{type}_{level}.mpa.txt"
	conda:
		"envs/bracken.yaml"
	log:
		log1 = LOG_DIR + "combine_result/{type}_{level}.combine_bracken.log",
		log2 = LOG_DIR + "combine_result/{type}_{level}.combine_mpa.log"
	message:
		"Executing Combine on the following files {input}."
	priority: 1
	shell:
		"""
		combine_bracken_outputs.py --files {params.input1} -o {output.op1} &> {log.log1}
		combine_mpa.py -i {params.input2} -o {output.op2} &> {log.log2}
		"""
