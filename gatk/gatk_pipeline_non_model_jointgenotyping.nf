#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/* * PIPELINE PARAMETERS 
 */
params.reads          = "/path/to/reads/*_{1,2}.fastq.gz" 
params.ref            = "/path/to/reference.fasta"
params.ref_index      = "/path/to/reference.fasta.*"       
params.interval_lists = "/path/to/polished_*.list" 
params.outdir         = "./results"
params.fastp_bin      = "fastp" 
params.cohort_name    = "cohort_joint_called" // Name for final file

// ------------------------------------------------------------------
// PROCESS 1: FASTP (Trimming) - Unchanged from single sample
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
// PROCESS 2: ALIGNMENT - Unchanged from single sample
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
    def java_mem = avail_mem - 8 
    if (java_mem < 4) java_mem = 4

    """
    gatk --java-options "-Xmx${java_mem}g" FastqToSam \
        --FASTQ $r1 --FASTQ2 $r2 \
        --OUTPUT unmapped.bam \
        --READ_GROUP_NAME ${sample_id} --SAMPLE_NAME ${sample_id} \
        --LIBRARY_NAME TruSeq --PLATFORM_UNIT NovaSeqX \
        --PLATFORM Illumina --SEQUENCING_CENTER DUGC

    set -o pipefail
    
    gatk --java-options "-Xmx${java_mem}g" SamToFastq \
        -I unmapped.bam -F /dev/stdout \
        --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 \
        --INTERLEAVE true --NON_PF true \
    | bwa-mem2 mem -M -t ${task.cpus} -p ${ref} /dev/stdin \
    | gatk --java-options "-Xmx${java_mem}g" MergeBamAlignment \
        --ALIGNED_BAM /dev/stdin --UNMAPPED_BAM unmapped.bam \
        --OUTPUT piped.bam -R ${ref} \
        --CREATE_INDEX true --ADD_MATE_CIGAR true \
        --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true \
        --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS

    gatk --java-options "-Xmx${java_mem}g" MarkDuplicatesWithMateCigar \
        --INPUT piped.bam \
        --OUTPUT ${sample_id}.md.bam \
        --METRICS_FILE ${sample_id}.metrics.txt \
        --MINIMUM_DISTANCE 300 --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --CREATE_INDEX true
    """
}

// ------------------------------------------------------------------
// PROCESS 3: HAPLOTYPE CALLER - Unchanged for sinsle sample
// ------------------------------------------------------------------
process HAPLOTYPE_CALLER {
    tag "${sample_id}-${interval_list.baseName}"
    publishDir "${params.outdir}/gvcfs_parts/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    each path(interval_list) 
    path ref
    path ref_dict
    path ref_idx

    output:
    tuple val(sample_id), path("${sample_id}.${interval_list.baseName}.g.vcf.gz"), emit: gvcfs

    script:
    def java_mem = task.memory.toGiga() - 2
    """
    gatk --java-options "-Xmx${java_mem}g" HaplotypeCaller \
        -R ${ref} -I ${bam} \
        -O ${sample_id}.${interval_list.baseName}.g.vcf.gz \
        -L ${interval_list} -ERC GVCF \
        --native-pair-hmm-threads ${task.cpus}
    """
}

// ------------------------------------------------------------------
// PROCESS 4: MERGE SINGLE SAMPLE GVCFS
// ------------------------------------------------------------------
// This fulfills your request to create single sample files before joint calling
process MERGE_SINGLE_SAMPLE_GVCFS {
    tag "$sample_id"
    publishDir "${params.outdir}/gvcfs_merged/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(gvcfs)
    path ref_dict

    output:
    // IMPORTANT: We output BOTH the vcf and the tbi (index). 
    // GenomicsDBImport needs the index to read these files randomly.
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: merged_gvcf

    script:
    """
    ls *.g.vcf.gz > input.list

    gatk MergeVcfs \
        -I input.list \
        -O ${sample_id}.g.vcf.gz \
        -D ${ref_dict}
    """
}

// ------------------------------------------------------------------
// PROCESS 5: GENOMICS DB IMPORT (Parallelized by Interval)
// ------------------------------------------------------------------
// We import ALL samples into a DB, but we do it separately for each interval list
process GENOMICS_DB {
    tag "${interval_list.baseName}"
    // No publishDir needed really, this is an intermediate folder structure
    
    input:
    path interval_list            // One interval file
    path gvcfs                    // List of ALL sample GVCFs
    path indices                  // List of ALL sample indices
    
    output:
    // Output is a directory (the database)
    tuple val(interval_list.baseName), path("gendb_${interval_list.baseName}"), emit: db_dir

    script:
    def java_mem = task.memory.toGiga() - 2
    // We construct a sample map or pass files via arguments. 
    // With many files, creating a sample map is safer than a long command line.
    """
    # Create sample map: sample_name <tab> path
    for f in *.g.vcf.gz; do
        # Extract sample name (assumed to be filename up to .g.vcf)
        name=\$(basename \$f .g.vcf.gz)
        echo -e "\$name\t\$f" >> sample_map.txt
    done

    gatk --java-options "-Xmx${java_mem}g" GenomicsDBImport \
        --genomicsdb-workspace-path gendb_${interval_list.baseName} \
        --sample-name-map sample_map.txt \
        -L ${interval_list} \
        --reader-threads ${task.cpus} \
        --batch-size 50 
    """
}

// ------------------------------------------------------------------
// PROCESS 6: JOINT GENOTYPING (Parallelized by Interval)
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
    gatk --java-options "-Xmx${java_mem}g" GenotypeGVCFs \
        -R ${ref} \
        -V gendb://${workspace_dir} \
        -O ${params.cohort_name}.${interval_name}.vcf.gz \
        -A AlleleFraction -A FisherStrand -A StrandOddsRatio -A ExcessHet
    """
}

// ------------------------------------------------------------------
// PROCESS 7: MERGE COHORT VCFS
// ------------------------------------------------------------------
// Gather the parallel chunks into one final file
process MERGE_COHORT_FINAL {
    publishDir "${params.outdir}/joint_genotyped", mode: 'copy'

    input:
    path vcfs // List of all chunked VCFs
    path ref_dict

    output:
    tuple path("${params.cohort_name}.vcf.gz"), path("${params.cohort_name}.vcf.gz.tbi"), emit: final_vcf

    script:
    """
    ls *.vcf.gz > input.list

    gatk MergeVcfs \
        -I input.list \
        -O ${params.cohort_name}.vcf.gz \
        -D ${ref_dict}
    """
}

// ------------------------------------------------------------------
// PROCESS 8: HARD FILTERING (On Cohort)
// ------------------------------------------------------------------
process FILTER_COHORT {
    publishDir "${params.outdir}/joint_genotyped", mode: 'copy'

    input:
    tuple path(vcf), path(idx)
    path ref
    path ref_dict
    path ref_idx

    output:
    tuple path("${params.cohort_name}.snps.filtered.vcf.gz"), path("${params.cohort_name}.indels.filtered.vcf.gz"), emit: filtered_final

    script:
    def base = params.cohort_name
    """
    # 1. Select SNPs
    gatk SelectVariants -R $ref -V $vcf -select-type SNP -O snps.vcf.gz
    
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
        --filter-name "ReadPosRankSum-8" --filter-expression "ReadPosRankSum < -8.0"

    # 3. Select Indels
    gatk SelectVariants -R $ref -V $vcf -select-type INDEL -O indels.vcf.gz

    # 4. Filter Indels
    gatk VariantFiltration \
        -R $ref -V indels.vcf.gz \
        -O ${base}.indels.filtered.vcf.gz \
        --filter-name "QD2"   --filter-expression "QD < 2.0" \
        --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
        --filter-name "FS200" --filter-expression "FS > 200.0" \
        --filter-name "ReadPosRankSum-20" --filter-expression "ReadPosRankSum < -20.0"
    """
}

// ------------------------------------------------------------------
// WORKFLOW
// ------------------------------------------------------------------
workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads)
    intervals_ch  = Channel.fromPath(params.interval_lists) 

    ref_ch = file(params.ref)
    ref_dict_ch = file(params.ref.replace(".fasta", ".dict").replace(".fa", ".dict"))
    ref_idx_ch = Channel.fromPath(params.ref_index).collect() 

    // --- UPSTREAM (Per Sample) ---
    FASTP(read_pairs_ch)
    ALIGN_AND_PREPROCESS(FASTP.out.trimmed_reads, ref_ch, ref_idx_ch, ref_dict_ch)
    HAPLOTYPE_CALLER(ALIGN_AND_PREPROCESS.out.bam_tuple, intervals_ch, ref_ch, ref_dict_ch, ref_idx_ch)
    MERGE_SINGLE_SAMPLE_GVCFS(HAPLOTYPE_CALLER.out.gvcfs.groupTuple(), ref_dict_ch)

    // --- THE GATHERING (Wait for ALL samples) ---
    // .collect() waits until every single sample is processed, 
    // then outputs a list [file1, file2, file3...]
    
    all_gvcfs_ch = MERGE_SINGLE_SAMPLE_GVCFS.out.merged_gvcf.map{ it[1] }.collect()
    all_tbis_ch  = MERGE_SINGLE_SAMPLE_GVCFS.out.merged_gvcf.map{ it[2] }.collect()

    // --- DOWNSTREAM (Joint) ---
    // We cross the intervals with the FULL list of gVCFs
    // Input: [Interval_File, [List_of_All_GVCFs], [List_of_All_TBIs]]
    
    GENOMICS_DB(intervals_ch, all_gvcfs_ch, all_tbis_ch)
    
    GENOTYPE_COHORT(GENOMICS_DB.out.db_dir, ref_ch, ref_dict_ch, ref_idx_ch)
    
    // Merge the scattered joint VCFs back into one final file
    MERGE_COHORT_FINAL(GENOTYPE_COHORT.out.vcf_parts.collect(), ref_dict_ch)
    
    FILTER_COHORT(MERGE_COHORT_FINAL.out.final_vcf, ref_ch, ref_dict_ch, ref_idx_ch)
}
