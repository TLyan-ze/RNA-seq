#################################################################################################################
###yanzeqin
###2020/10/29
#################################################################################################################
##根据input的数据量，来控制该rule的线程数
##隔段时间统计盘阵资源，如果盘阵资源较少则把所有任务挂起，直到盘阵空间足够再恢复
#################################################################################################################

import snakemake

configfile: "config.yaml"
print (config['samples'])




rule all:
    input:
        expand("/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}.merged.fastq",
               sample=config["samples"])


#########################################FQ######################################################################
rule FQ:
	input:
		"_1.fq.gz","_2.fq.gz"
	output:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/01.fastp/{sample}_1.fq.gz",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/01.fastp/{sample}_2.fq.gz",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/01.fastp/{sample}.report.html",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/01.fastp/{sample}.report.json"
	shell:
		"fastp -w 8 -i {input[0]} -I {input[1]}  \
        	--adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        	--adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        	--detect_adapter_for_pe -q 20 -u 20 \
        	--length_required 50 -n 2 -y -c -p \
        	--disable_trim_poly_g \
        	-o {output[0]} -O {output[1]} \
        	-h {output[2]} -j {output[3]}"



rule SOAPnuke_filter:
    input:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/01.fastp/{sample}_1.fq.gz",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/01.fastp/{sample}_2.fq.gz"
    output:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/02.SOA/{sample}_1.fq.gz",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/02.SOA/{sample}_2.fq.gz"
    params:
        rg="/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/02.SOA/",
        cg="{sample}_1.fg.gz",
        dg="{sample}_2.fg.gz"
    shell:
        "/ldfssz1/ST_INFECTION/F16ZQSB1SY3030_Salmonella/USER/lindechun/bin/SOAPnuke filter -T 8 \
        -1 {input[0]} \
        -2 {input[1]} \
        -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        -E 50 -l 20 -q 0.2 -n 0.02 -d -Q 2 -G \
        -C {params.cg} -D {params.dg} \
        -o {params.rg}"


rule prinseq:
    input:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/02.SOA/{sample}_1.fq.gz",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/02.SOA/{sample}_2.fq.gz"
    output:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/03.SOA/{sample}",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/03.SOA/{sample}_good_out_R1.fastq",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/03.SOA/{sample}_good_out_R2.fastq"
    shell:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/prinseq++_1.2/bin/prinseq++ \
        -threads 8 \
        -fastq {input[0]} \
        -fastq2 {input[1]} \
        -lc_entropy=0.5 -lc_dust=0.5 \
        -out_name {output[0]} "


#############################################sortmeRNA#################################################
rule rmrRNA_1:
    input:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/03.SOA/{sample}_good_out_R1.fastq",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/01.QC/03.SOA/{sample}_good_out_R2.fastq"
    output:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}.merged.fastq"
    shell:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/scripts/merge-paired-reads.sh {input[0]} {input[1]} {output}"


rule rmrRNA_2:
    input:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}.merged.fastq"
    output:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}_rmRNA.fastq"
    params:
        r1="{sample}",
        r2="{sample}_rmRNA"
    shell:
        "sortmerna -a 5 -m 10240 \
        --ref $rRNAdb/rfam-5.8s-database-id98.fasta,\
        $rRNAindex/rfam-5.8s-db:$rRNAdb/rfam-5s-database-id98.fasta,\
        $rRNAindex/rfam-5s-db:$rRNAdb/silva-arc-16s-id95.fasta,\
        $rRNAindex/silva-arc-16s-db:$rRNAdb/silva-arc-23s-id98.fasta,\
        $rRNAindex/silva-arc-23s-db:$rRNAdb/silva-bac-16s-id90.fasta,\
        $rRNAindex/silva-bac-16s-db:$rRNAdb/silva-bac-23s-id98.fasta,\
        $rRNAindex/silva-bac-23s-db:$rRNAdb/silva-euk-18s-id95.fasta,\
        $rRNAindex/silva-euk-18s-db:$rRNAdb/silva-euk-28s-id98.fasta,\
        $rRNAindex/silva-euk-28s

        /hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/sortmerna \
        -a 5 -m 10240 \
        --ref /hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,\
        /hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/rfam-5.8s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,\
        /hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/rfam-5s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,\
        /hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-arc-16s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-arc-23s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-bac-16s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-bac-23s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-euk-18s-db:/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/index/silva-euk-28s
        --reads {input} \
        --aligned {params.r1} \
        --other {params.r2} \
        --num_alignments 1 --fastx --paired_out --log -v"
    

rule rrmrRNA_3:
    input:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}_rmRNA.fastq"
    output:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}_1_rmRNA.fastq",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}_2_rmRNA.fastq"
    shell:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/share/sortmerna-2.1b/scripts/unmerge-paired-reads.sh \
        $outdir/$sample/$batch/02.rmRNA/${sample}_rmRNA.fastq \
        $outdir/$sample/$batch/02.rmRNA/${sample}_1_rmRNA.fastq \
        $outdir/$sample/$batch/02.rmRNA/${sample}_2_rmRNA.fastq"


####################################kraken2###################################################
rule kraken2:
    input:
        ""
    output:
    shell:
        "kraken2 --db $K2DB --threads 5 \
        --report $outdir/$sample/Merged/03.Kraken2/${sample}.kraken2.report \
        --paired $outdir/$sample/Merged/02.rmRNA/${sample}_1_rmRNA.fastq \
        $outdir/$sample/Merged/02.rmRNA/${sample}_2_rmRNA.fastq \
        > $outdir/$sample/Merged/03.Kraken2/${sample}.kraken2"


 


rule braken:
    input:
    output:
    shell:
        "bracken -d $K2DB -i $outdir/$sample/Merged/03.Kraken2/${sample}.kraken2.report \
        -o $outdir/$sample/Merged/03.Kraken2/Bracken/${sample}.Family.bracken \
        -w $outdir/$sample/Merged/03.Kraken2/Bracken/${sample}.Family.bracken.report -l F |"
        "bracken -d $K2DB -i $outdir/$sample/Merged/03.Kraken2/${sample}.kraken2.report \
        -o $outdir/$sample/Merged/03.Kraken2/Bracken/${sample}.Genus.bracken \
        -w $outdir/$sample/Merged/03.Kraken2/Bracken/${sample}.Genus.bracken.report -l G |"
        "bracken -d $K2DB -i $outdir/$sample/Merged/03.Kraken2/${sample}.kraken2.report \
        -o $outdir/$sample/Merged/03.Kraken2/Bracken/${sample}.Species.bracken \
        -w $outdir/$sample/Merged/03.Kraken2/Bracken/${sample}.Species.bracken.report -l S"




###########################################组装&&&&####################################################################
rule trinity:
    input:
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}_1_rmRNA.fastq",
        "/hwfssz5/ST_INFECTION/GlobalDatabase/user/yanzeqin/ggg/{sample}/02.rmRNA/{sample}_2_rmRNA.fastq"
    output:
        ""
    shell:
        "Trinity --seqType fq --left $outdir/$sample/Merged/02.rmRNA/${sample}_1_rmRNA.fastq  \
        --right $outdir/$sample/Merged/02.rmRNA/${sample}_2_rmRNA.fastq  \
        --min_contig_length 1000  --CPU 10 --max_memory 200G --no_version_check \
         --output $outdir/$sample/Merged/04.Trinity/"




rule cd_hit-est:
    input:
    output:
    shell:
        "$cdhitest -i $outdir/$sample/$batch/results/04.Trinity/Trinity.fasta \
        -o $outdir/$sample/$batch/results/05.rmDuplications/${sample}.unique.contigs.fasta \
        -M 10240 -T 5 -c 0.95 -s 0.9 -aS 0.9 -g 1 &&"

    







##############################################mapping back&&&coverage#######################################
rule bowtie:
    input:
    output:

    shell:
        ""





rule so