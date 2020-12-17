#####################################################
###yanzeqin
###2020/10/29
#######################################################
# read config info into this namespace




configfile: "config.yaml"


SAMPLES = {}
with open(config['params']['samples'], 'rb') as f:
    for line in f:
        parts = line.decode('utf-8').split()
        key = parts[0]
        SAMPLES[key] = parts[1]
        print(key)
print(SAMPLES)





rule all:
	input:
		expand("result/{sample}/04.metaSPAdes/scaffolds.fasta",sample=SAMPLES),
		expand("result/{sample}/03.Kraken2/Bracken/{sample}.Family.bracken",sample=SAMPLES)



##############################################01.QC############################################
rule FQ:
	input:
		lambda wildcards: f"{SAMPLES[wildcards.sample]}1.fq.gz",
		lambda wildcards: f"{SAMPLES[wildcards.sample]}2.fq.gz"
	output:
		"result/{sample}/01.QC/0101.fastp/{sample}_1.fq.gz",
		"result/{sample}/01.QC/0101.fastp/{sample}_2.fq.gz",
		"result/{sample}/01.QC/0101.fastp/{sample}.report.html",
		"result/{sample}/01.QC/0101.fastp/{sample}.report.json",
		"result/{sample}/raw.stats"
	shell:
		"""
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/softwares1/ALL_tools/Fastp/fastp -w 8 -i {input[0]} -I {input[1]} \
		--adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
		--adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
		--detect_adapter_for_pe -q 20 -u 20 \
		--length_required 50 -n 2 -y -c -p \
		--disable_trim_poly_g \
		-o {output[0]} -O {output[1]} \
		-h {output[2]} -j {output[3]}
		seqkit stats -j 2 -a  {input[0]}  > {output[4]}
		"""



##SOAPnuke filter
rule SOAPnuke_filter:
	input:
		"result/{sample}/01.QC/0101.fastp/{sample}_1.fq.gz",
		"result/{sample}/01.QC/0101.fastp/{sample}_2.fq.gz"
	output:
		"result/{sample}/01.QC/0102.soapnuke/{sample}_1.fq.gz",
		"result/{sample}/01.QC/0102.soapnuke/{sample}_2.fq.gz",
		"result/{sample}/fastp.stats"
	params:
		rg="result/{sample}/01.QC/0102.soapnuke/",
		cg="{sample}_1.fq.gz",
		dg="{sample}_2.fq.gz"
	shell:
		"""
		/ldfssz1/ST_INFECTION/F16ZQSB1SY3030_Salmonella/USER/lindechun/bin/SOAPnuke filter -T 8 \
		-1 {input[0]} -2 {input[1]} \
		-f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
		-r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
		-E 50 -l 20 -q 0.2 -n 0.02 -d -Q 2 -G \
		-C {params.cg} -D {params.dg} \
		-o {params.rg}
		seqkit stats -j 2 -a  {input[0]} > {output[2]}
		"""




rule prinseq:
	input:
		"result/{sample}/01.QC/0102.soapnuke/{sample}_1.fq.gz",
		"result/{sample}/01.QC/0102.soapnuke/{sample}_2.fq.gz"
	output:
		"result/{sample}/01.QC/0103.prinseq/{sample}_good_out_R1.fastq",
		"result/{sample}/01.QC/0103.prinseq/{sample}_good_out_R2.fastq",
		"result/{sample}/soapnuke.stats"
	params:
		"result/{sample}/01.QC/0103.prinseq/{sample}"
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/prinseq++_1.2/bin/prinseq++ -threads 8 -fastq {input[0]} -fastq2 {input[1]} -lc_entropy=0.5 -lc_dust=0.5 -out_name {params[0]}
		seqkit stats -j 2 -a  {input[0]} > {output[2]}
		'''



rule rmrRNA_1:
	input:
		"result/{sample}/01.QC/0103.prinseq/{sample}_good_out_R1.fastq",
		"result/{sample}/01.QC/0103.prinseq/{sample}_good_out_R2.fastq"
	output:
		"result/{sample}/01.QC/0104.sortmeRNA/{sample}.merged.fastq",
		"result/{sample}/prinseq.stats"
	shell:
		"""
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/scripts/merge-paired-reads.sh {input[0]} {input[1]} {output[0]}
		seqkit stats -j 2 -a  {input[0]} > {output[1]}
		"""




rule rmrRNA_2:
	input:
		"result/{sample}/01.QC/0104.sortmeRNA/{sample}.merged.fastq"
	output:
		"result/{sample}/02.rmRNA/{sample}_rmRNA.fastq"
	params:
		rmRNA0="result/{sample}/02.rmRNA/{sample}",
		rmRNA1="result/{sample}/02.rmRNA/{sample}_rmRNA"
	shell:
		"/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/sortmerna -a 5 -m 10240 \
		--ref /hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/rfam-5.8s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/rfam-5s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-arc-16s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-arc-23s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-bac-16s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-bac-23s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-euk-18s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-euk-28s \
		--reads {input} --aligned {params.rmRNA0} --other {params.rmRNA1} --num_alignments 1 --fastx --paired_out --log -v"


rule rmrRNA_3:
	input:
		"result/{sample}/02.rmRNA/{sample}_rmRNA.fastq"
	output:
		"result/{sample}/02.rmRNA/{sample}_1_rmRNA.fastq",
		"result/{sample}/02.rmRNA/{sample}_2_rmRNA.fastq"
	shell:
		"""
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/scripts/unmerge-paired-reads.sh {input} {output[0]} {output[1]}
		"""



#################################################kraken2###################################################################
rule kraken2:
	input:
		"result/{sample}/02.rmRNA/{sample}_1_rmRNA.fastq",
		"result/{sample}/02.rmRNA/{sample}_2_rmRNA.fastq"
	output:
		"result/{sample}/03.kraken2/{sample}.kraken2.report",
		"result/{sample}/03.kraken2/{sample}.kraken2",
		"result/{sample}/rmrRNA.stats"
	shell:
		"""
		kraken2 --db /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/K2DB/ncbik2db \
		--threads 5 --report {output[0]}  \
		--paired {input[0]} {input[1]} > {output[1]}
		seqkit stats -j 2 -a  {input[0]} > {output[2]}
		"""
		



rule barcken:
	input:
		"result/{sample}/03.kraken2/{sample}.kraken2.report"
	output:
		"result/{sample}/03.Kraken2/Bracken/{sample}.Family.bracken",
		"result/{sample}/03.Kraken2/Bracken/{sample}.Genus.bracken",
		"result/{sample}/03.Kraken2/Bracken/{sample}.Species.bracken",
		"result/{sample}/03.Kraken2/Bracken/{sample}.Family.bracken.report",
                "result/{sample}/03.Kraken2/Bracken/{sample}.Genus.bracken.report",
		"resul



t/{sample}/03.Kraken2/Bracken/{sample}.Species.bracken.report"

	shell:
		'''
		bracken -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/K2DB/ncbik2db -i {input} -o {output[0]} -w {output[3]} -l F
		bracken -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/K2DB/ncbik2db -i {input} -o {output[1]} -w {output[4]} -l G
		bracken -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/K2DB/ncbik2db -i {input} -o {output[2]} -w {output[5]} -l S 
		'''

##############################################Trinity##assemble#########################################################################
######trinity


rule spades:
	input:
		"result/{sample}/02.rmRNA/{sample}_1_rmRNA.fastq",
		"result/{sample}/02.rmRNA/{sample}_2_rmRNA.fastq"
	output:
		"result/{sample}/04.metaSPAdes/scaffolds.fasta"
	params:
		"result/{sample}/04.metaSPAdes"
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/SPAdes-3.14.1/bin/metaspades.py \
		-1 {input[0]} \
		-2 {input[1]} \
		-t 10 -m 300 -o {params[0]}
		'''






###########################################聚类###初筛###注释#################################################

rule cdhitest:
	input:
		"result/{sample}/04.metaSPAdes/scaffolds.fasta"
	output:
		"result/{sample}/05.rmDup/{sample}.unique.contigs.fasta"
	shell:
		"/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/cd-hit-4.8.1/bin/cd-hit-est -i {input} -o {output} -M 10240 -T 5 -c 0.95 -s 0.9 -aS 0.9 -g 1 "




rule blastx:
	input:
		"result/{sample}/05.rmDup/{sample}.unique.contigs.fasta"
	output:
		"result/{sample}/05.rmDup/{sample}.blast"
	shell:
		"diamond blastx -d -q {input} -k 1 -f 6 -e 0.00001 -o {output} "


rule extract_id:
	input:
		"result/{sample}/05.rmDup/{sample}.blast"
	output:
		"result/{sample}/05.rmDup/{sample}.id.txt"
	shell:
		" cat {input} | awk '{print $1}' > {output}"



rule id_extract_fastq:
	input:
		"result/{sample}/05.rmDup/{sample}.unique.contigs.fasta",
		"result/{sample}/05.rmDup/{sample}.id.txt"
	output:
		"result/{sampel}/05.rmDup/{sample}.blast.fasta"
	shell:
		"""
		seqkit grep -f {input[1]} {input[0]} > {output[0]}
		"""


rule CAT_1:
	input:
		"result/{sample}/05.rmDup/{sample}.unique.contigs.fasta"
	output:
		"result/{sample}/06.CAT/{sample}.contig2classification.txt",
		"result/{sample}/06.CAT/{sample}.ORF2LCA.txt"
	params:
		"result/{sample}/06.CAT/{sample}"
	threads:4
 	shell:
                '''
                /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT contigs -n 4 --force --block_size 10.0 --index_chunks 1 -c {input} -o {params[0]} --no_stars -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/nr_database/CAT_nr/ -t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --path_to_prodigal /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Prodigal-2.6.3/prodigal --path_to_diamond /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/diamond
                '''


rule CAT_2:
	input:
		"result/{sample}/06.CAT/{sample}.contig2classification.txt",
		"result/{sample}/06.CAT/{sample}.ORF2LCA.txt"
	output:
		"result/{sample}/06.CAT/{sample}.contig.tax",
		"result/{sample}/06.CAT/{sample}.contig.official.tax",
		"result/{sample}/06.CAT/{sample}.ORF.tax"
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT add_names -i {input[0]} -o {output[1]} -t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --only_official
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT add_names -i {input[0]} -o {output[0]} -t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT add_names -i {input[1]} -o {output[2]} -t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --only_official
		'''

####
rule CAT_3:
	input:
		"result/{sample}/06.CAT/{sample}.contig.official.tax",
		"result/{sample}/05.rmDup/{sample}.unique.contigs.fasta"
	output:
		"result/{sample}/06.CAT/{sample}.summary.txt"
	shell:
		"/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT summarise -c {input[1]} -i {input[0]} -o {output}"




########################################mapping back&&&coverage#####################################
rule index_MappingBack:
	input:
		"result/{sample}/04.metaSPAdes/scaffolds.fasta"
	output:
		"result/{sample}/04.metaSPAdes/{sample}_spades.1.bt2"
	params:
		index1="result/{sample}/04.metaSPAdes/{sample}_spades"
	shell:
		"bowtie2-build {input}  {params.index1}"



rule mapping_MappingBack:
	input:
		"result/{sample}/04.metaSPAdes/{sample}_spades.1.bt2",
		"result/{sample}/02.rmRNA/{sample}_1_rmRNA.fastq",
		"result/{sample}/02.rmRNA/{sample}_2_rmRNA.fastq"
	output:
		"result/{sample}/07.sam/{sample}.sam"
	params:
		index2="result/{sample}/04.metaSPAdes/{sample}_spades"
	shell:
		"bowtie2 -p 8 -x {params.index2} -1 {input[1]} -2 {input[2]} -S {output}"



rule soapcoverage_mappingBack:
	input:
		"result/{sample}/07.sam/{sample}.sam",
		"result/{sample}/04.metaSPAdes/scaffolds.fasta"
	output:
		"result/{sample}/08.coverage/{sample}_sofa.txt"
	shell:
		"/ldfssz1/ST_INFECTION/P18Z10200N0164_Resistance_PN/User/zhaohailong/software/soap.coverage \
		-cvg -sam -refsingle {input[1]} -i {input[0]} -o {output} "



rule stats_all:
	input:
		raw=expand("result/{sample}/raw.stats",sample=SAMPLES),
		fastp=expand("result/{sample}/fastp.stats",sample=SAMPLES),
		soapnuke=expand("result/{sample}/soapnuke.stats",sample=SAMPLES),
		prinseq=expand("result/{sample}/prinseq.stats",sample=SAMPLES),
		rmrRNA=expand("result/{sample}/rmrRNA.stats",sample=SAMPLES)
	outut:
	shell:

