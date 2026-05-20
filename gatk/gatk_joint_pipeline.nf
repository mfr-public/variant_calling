#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/* * PIPELINE PARAMETERS 
 */
params.reads          = "/scratch/pawsey0149/mrichardson/raw/fastqs/compiled/*_{1,2}.fastq.gz" 
params.ref            = "/scratch/pawsey0149/mrichardson/genome/polished.fasta"
params.ref_index      = "/scratch/pawsey0149/mrichardson/genome/polished.fasta.*"       
params.interval_lists = "/scratch/pawsey0149/mrichardson/genome/intervals/*.interval_list" 
params.outdir         = "./results"
params.tmp_dir        = "/scratch/pawsey0149/mrichardson/ng_all/gatk_tmp" 
params.fastp_bin      = "fastp" 
params.cohort_name    = "cohort_joint_called" 

// ------------------------------------------------------------------
// PROCESS 1: FASTP 
// ------------------------------------------------------------------
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_trim.fastq.gz"), path("${sample_id}_2_trim.fastq.gz"), emit: trimmed_reads
    path "${sample_id}.fastp.json"

    script:
    """
    ${params.fastp_bin} \
        -i ${reads[0]} -I ${reads[1]} \
        -o ${sample_id}_1_trim.fastq.gz -O ${sample_id}_2_trim.fastq.gz \
        --detect_adapter_for_pe \
        --trim_tail1 1 --trim_tail2 1 \
        -g --poly_g_min_len 1 \
        -q 20 -u 90 -l 36 \
        -j ${sample_id}.fastp.json
    """
}

// ------------------------------------------------------------------
// PROCESS 2: ALIGNMENT & PRE-PROCESSING
// ------------------------------------------------------------------
process ALIGN_AND_PREPROCESS {
    tag "$sample_id"
    publishDir "${params.outdir}/bams/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path ref
    path ref_idx 
    path ref_dict 

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), emit: bam_tuple

    script:
    def avail_mem = task.memory.toGiga()
    def mem_single = avail_mem - 8
    if (mem_single < 8) mem_single = 8
    def mem_pipe = 4 // HARDCAP for streaming tools

    """
    mkdir -p ${params.tmp_dir}

    gatk --java-options "-Xmx${mem_single}g" FastqToSam \
        --FASTQ $r1 --FASTQ2 $r2 \
        --OUTPUT unmapped.bam \
        --READ_GROUP_NAME ${sample_id} --SAMPLE_NAME ${sample_id} \
        --LIBRARY_NAME TruSeq --PLATFORM_UNIT NovaSeqX \
        --PLATFORM Illumina --SEQUENCING_CENTER NG \
        --TMP_DIR ${params.tmp_dir}

    set -o pipefail
    
    gatk --java-options "-Xmx${mem_pipe}g" SamToFastq \
        -I unmapped.bam -F /dev/stdout \
        --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 \
        --INTERLEAVE true --NON_PF true \
        --TMP_DIR ${params.tmp_dir} \
    | bwa-mem2 mem -M -t ${task.cpus} -p ${ref} /dev/stdin \
    | gatk --java-options "-Xmx${mem_pipe}g" MergeBamAlignment \
        --ALIGNED_BAM /dev/stdin --UNMAPPED_BAM unmapped.bam \
        --OUTPUT piped.bam -R ${ref} \
        --CREATE_INDEX true --ADD_MATE_CIGAR true \
        --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true \
        --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS \
        --TMP_DIR ${params.tmp_dir}

    gatk --java-options "-Xmx${mem_single}g" MarkDuplicatesWithMateCigar \
        --INPUT piped.bam \
        --OUTPUT ${sample_id}.md.bam \
        --METRICS_FILE ${sample_id}.metrics.txt \
        --MINIMUM_DISTANCE 300 --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --CREATE_INDEX true \
        --TMP_DIR ${params.tmp_dir}
    """
}

// ------------------------------------------------------------------
// PROCESS 3: HAPLOTYPE CALLER
// ------------------------------------------------------------------
process HAPLOTYPE_CALLER {
    tag "${sample_id}-${interval_list.baseName}"
    // Optional: Turn off publishDir here if you dont want 6,600 files in your results folder
    publishDir "${params.outdir}/gvcfs_parts/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    each path(interval_list) 
    path ref
    path ref_dict
    path ref_idx

    output:
    // EXPLICITLY outputting the interval_name and index to enable Scatter-to-Scatter bypassing
    tuple val(sample_id), val(interval_list.baseName), path("${sample_id}.${interval_list.baseName}.g.vcf.gz"), path("${sample_id}.${interval_list.baseName}.g.vcf.gz.tbi"), emit: gvcfs

    script:
    def java_mem = task.memory.toGiga() - 2
    """
    mkdir -p ${params.tmp_dir}

    gatk --java-options "-Xmx${java_mem}g" HaplotypeCaller \
        -R ${ref} -I ${bam} \
        -O ${sample_id}.${interval_list.baseName}.g.vcf.gz \
        -L ${interval_list} -ERC GVCF \
        --native-pair-hmm-threads ${task.cpus} \
        --tmp-dir ${params.tmp_dir}
    """
}

// ------------------------------------------------------------------
// PROCESS 4: MERGE SINGLE SAMPLE GVCFS (Archive Track)
// ------------------------------------------------------------------
process MERGE_SINGLE_SAMPLE_GVCFS {
    tag "$sample_id"
    publishDir "${params.outdir}/gvcfs_merged/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(gvcfs)
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: merged_gvcf

    script:
    def input_args = gvcfs.collect { "-I ${it}" }.join(' ') 
    """
    mkdir -p ${params.tmp_dir}

    gatk MergeVcfs \
        ${input_args} \
        -O ${sample_id}.g.vcf.gz \
        -D ${ref_dict} \
        --tmp-dir ${params.tmp_dir}
    """
}

// ------------------------------------------------------------------
// PROCESS 5: GENOMICS DB IMPORT (Fast Track)
// ------------------------------------------------------------------
process GENOMICS_DB {
    tag "${interval_name}"
    
    input:
    tuple val(interval_name), path(gvcfs), path(indices), path(interval_list)
    
    output:
    tuple val(interval_name), path("gendb_${interval_name}"), emit: db_dir

    script:
    def java_mem = task.memory.toGiga() - 2
    """
    mkdir -p ${params.tmp_dir}

    # Dynamically generate the sample map locally using the scattered chunks
    for f in *.g.vcf.gz; do
        # Extract true sample name by trimming the interval and extension
        name=\$(basename \$f .${interval_name}.g.vcf.gz)
        echo -e "\$name\t\$f" >> sample_map.txt
    done

    gatk --java-options "-Xmx${java_mem}g" GenomicsDBImport \
        --genomicsdb-workspace-path gendb_${interval_name} \
        --sample-name-map sample_map.txt \
        -L ${interval_list} \
        --reader-threads ${task.cpus} \
        --batch-size 50 \
        --tmp-dir ${params.tmp_dir}
    """
}

// ------------------------------------------------------------------
// PROCESS 6: JOINT GENOTYPING
// ------------------------------------------------------------------
process GENOTYPE_COHORT {
    tag "${interval_name}"
    
    input:
    tuple val(interval_name), path(workspace_dir)
    path ref
    path ref_dict
    path ref_idx

    output:
    path "${params.cohort_name}.${interval_name}.vcf.gz", emit: vcf_parts

    script:
    def java_mem = task.memory.toGiga() - 2
    """
    mkdir -p ${params.tmp_dir}

    gatk --java-options "-Xmx${java_mem}g" GenotypeGVCFs \
        -R ${ref} \
        -V gendb://${workspace_dir} \
        -O ${params.cohort_name}.${interval_name}.vcf.gz \
        -A AlleleFraction -A FisherStrand -A StrandOddsRatio -A ExcessHet -A InbreedingCoeff \
        --tmp-dir ${params.tmp_dir}
    """
}

// ------------------------------------------------------------------
// PROCESS 7: MERGE COHORT VCFS
// ------------------------------------------------------------------
process MERGE_COHORT_FINAL {
    publishDir "${params.outdir}/joint_genotyped", mode: 'copy'

    input:
    path vcfs 
    path ref_dict

    output:
    tuple path("${params.cohort_name}.vcf.gz"), path("${params.cohort_name}.vcf.gz.tbi"), emit: final_vcf

    script:
    def input_args = vcfs.collect { "-I ${it}" }.join(' ')
    """
    mkdir -p ${params.tmp_dir}

    gatk MergeVcfs \
        ${input_args} \
        -O ${params.cohort_name}.vcf.gz \
        -D ${ref_dict} \
        --tmp-dir ${params.tmp_dir}
    """
}

// ------------------------------------------------------------------
// PROCESS 8: HARD FILTERING
// ------------------------------------------------------------------
process FILTER_COHORT {
    publishDir "${params.outdir}/joint_genotyped", mode: 'copy'

    input:
    tuple path(vcf), path(idx)
    path ref
    path ref_dict
    path ref_idx

    output:
    tuple path("${params.cohort_name}.snps.filtered.vcf.gz"), path("${params.cohort_name}.snps.filtered.vcf.gz.tbi"), \
          path("${params.cohort_name}.indels.filtered.vcf.gz"), path("${params.cohort_name}.indels.filtered.vcf.gz.tbi"), emit: filtered_final

    script:
    def base = params.cohort_name
    """
    mkdir -p ${params.tmp_dir}

    # 0. Pre-filter Excess Heterozygosity 
    gatk VariantFiltration \
        -R $ref -V $vcf \
        -O excesshet_filtered.vcf.gz \
        --filter-expression "ExcessHet > 54.69" \
        --filter-name "ExcessHet" \
        --tmp-dir ${params.tmp_dir}

    # 1. Select SNPs
    gatk SelectVariants -R $ref -V excesshet_filtered.vcf.gz -select-type SNP -O snps.vcf.gz --tmp-dir ${params.tmp_dir}
    
    # 2. Filter SNPs
    gatk VariantFiltration \
        -R $ref -V snps.vcf.gz \
        -O ${base}.snps.filtered.vcf.gz \
        --filter-name "QD2"             --filter-expression "QD < 2.0" \
        --filter-name "QUAL30"          --filter-expression "QUAL < 30.0" \
        --filter-name "SOR3"            --filter-expression "SOR > 3.0" \
        --filter-name "FS60"            --filter-expression "FS > 60.0" \
        --filter-name "MQ40"            --filter-expression "MQ < 40.0" \
        --filter-name "MQRankSum-12.5"  --filter-expression "MQRankSum < -12.5" \
        --filter-name "ReadPosRankSum-8" --filter-expression "ReadPosRankSum < -8.0" \
        --tmp-dir ${params.tmp_dir}

    # 3. Select Indels
    gatk SelectVariants -R $ref -V excesshet_filtered.vcf.gz -select-type INDEL -O indels.vcf.gz --tmp-dir ${params.tmp_dir}

    # 4. Filter Indels
    gatk VariantFiltration \
        -R $ref -V indels.vcf.gz \
        -O ${base}.indels.filtered.vcf.gz \
        --filter-name "QD2"   --filter-expression "QD < 2.0" \
        --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
        --filter-name "FS200" --filter-expression "FS > 200.0" \
        --filter-name "ReadPosRankSum-20" --filter-expression "ReadPosRankSum < -20.0" \
        --tmp-dir ${params.tmp_dir}
    """
}

// ------------------------------------------------------------------
// WORKFLOW LOGIC
// ------------------------------------------------------------------
workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads)
    intervals_ch  = Channel.fromPath(params.interval_lists) 

    ref_ch = file(params.ref)
    ref_dict_ch = file(params.ref.replace(".fasta", ".dict").replace(".fa", ".dict"))
    ref_idx_ch = Channel.fromPath(params.ref_index).collect() 

    // --- UPSTREAM ---
    FASTP(read_pairs_ch)
    ALIGN_AND_PREPROCESS(FASTP.out.trimmed_reads, ref_ch, ref_idx_ch, ref_dict_ch)
    HAPLOTYPE_CALLER(ALIGN_AND_PREPROCESS.out.bam_tuple, intervals_ch, ref_ch, ref_dict_ch, ref_idx_ch)

    // --- TRACK 1: ARCHIVE (Merge by Sample ID) ---
    // Shape: [sample_id, gvcf]
    hc_for_merge = HAPLOTYPE_CALLER.out.gvcfs
        .map { sample_id, interval_name, gvcf, tbi -> tuple(sample_id, gvcf) }
        .groupTuple()
    MERGE_SINGLE_SAMPLE_GVCFS(hc_for_merge, ref_dict_ch)

    // --- TRACK 2: FAST TRACK JOINT CALLING (Merge by Interval) ---
    // Shape: [interval_name, gvcf, tbi]
    hc_for_gendb = HAPLOTYPE_CALLER.out.gvcfs
        .map { sample_id, interval_name, gvcf, tbi -> tuple(interval_name, gvcf, tbi) }
        .groupTuple()
    
    // Map the raw interval files so we can join them natively
    intervals_mapped = intervals_ch.map { it -> tuple(it.baseName, it) }
    
    // Join HC outputs with the specific interval file
    // Input format to GENOMICS_DB: [interval_name, [List_of_GVCFs], [List_of_TBIs], interval_file]
    gendb_input = hc_for_gendb.join(intervals_mapped)

    GENOMICS_DB(gendb_input)
    GENOTYPE_COHORT(GENOMICS_DB.out.db_dir, ref_ch, ref_dict_ch, ref_idx_ch)
    MERGE_COHORT_FINAL(GENOTYPE_COHORT.out.vcf_parts.collect(), ref_dict_ch)
    FILTER_COHORT(MERGE_COHORT_FINAL.out.final_vcf, ref_ch, ref_dict_ch, ref_idx_ch)
}