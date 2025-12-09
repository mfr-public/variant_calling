#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/* * PIPELINE PARAMETERS 
 */
params.reads          = "/path/to/reads/*_{1,2}.fastq.gz" 
params.ref            = "/path/to/reference.fasta"
params.ref_index      = "/path/to/reference.fasta.*"       
params.interval_lists = "/path/to/polished_*.list" 
params.outdir         = "./results"
// If you want to use your specific binary, pass this param. Otherwise leave 'fastp'
params.fastp_bin      = "fastp" 

// ------------------------------------------------------------------
// PROCESS 1: FASTP (Trimming)
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
    // Matches logic from _FqTrm2VCFcluster
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
    def java_mem = avail_mem - 8 
    if (java_mem < 4) java_mem = 4

    """
    # 1. Convert Fastq to uBAM (Unmapped BAM)
    # Note: LIBRARY_NAME and PLATFORM hardcoded to match your script
    gatk --java-options "-Xmx${java_mem}g" FastqToSam \
        --FASTQ $r1 --FASTQ2 $r2 \
        --OUTPUT unmapped.bam \
        --READ_GROUP_NAME ${sample_id} \
        --SAMPLE_NAME ${sample_id} \
        --LIBRARY_NAME TruSeq \
        --PLATFORM_UNIT NovaSeqX \
        --PLATFORM Illumina \
        --SEQUENCING_CENTER DUGC

    # 2. Pipe: uBAM -> SamToFastq -> BWA-MEM2 -> MergeBamAlignment
    set -o pipefail
    
    gatk --java-options "-Xmx${java_mem}g" SamToFastq \
        -I unmapped.bam -F /dev/stdout \
        --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 \
        --INTERLEAVE true --NON_PF true \
    | bwa-mem2 mem -M -t ${task.cpus} -p ${ref} /dev/stdin \
    | gatk --java-options "-Xmx${java_mem}g" MergeBamAlignment \
        --ALIGNED_BAM /dev/stdin \
        --UNMAPPED_BAM unmapped.bam \
        --OUTPUT piped.bam \
        -R ${ref} \
        --CREATE_INDEX true --ADD_MATE_CIGAR true \
        --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true \
        --INCLUDE_SECONDARY_ALIGNMENTS true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ATTRIBUTES_TO_RETAIN XS

    # 3. MarkDuplicatesWithMateCigar
    gatk --java-options "-Xmx${java_mem}g" MarkDuplicatesWithMateCigar \
        --INPUT piped.bam \
        --OUTPUT ${sample_id}.md.bam \
        --METRICS_FILE ${sample_id}.metrics.txt \
        --MINIMUM_DISTANCE 300 \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --CREATE_INDEX true
    """
}

// ------------------------------------------------------------------
// PROCESS 3: HAPLOTYPE CALLER (SCATTERED)
// ------------------------------------------------------------------
process HAPLOTYPE_CALLER {
    tag "${sample_id}-${interval_list.baseName}"
    publishDir "${params.outdir}/gvcfs/${sample_id}", mode: 'copy'
    
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
        -R ${ref} \
        -I ${bam} \
        -O ${sample_id}.${interval_list.baseName}.g.vcf.gz \
        -L ${interval_list} \
        -ERC GVCF \
        --native-pair-hmm-threads ${task.cpus}
    """
}

// ------------------------------------------------------------------
// PROCESS 4: MERGE GVCFS
// ------------------------------------------------------------------
process MERGE_GVCFS {
    tag "$sample_id"
    publishDir "${params.outdir}/gvcfs/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(gvcfs)
    path ref_dict

    output:
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
// PROCESS 5: GENOTYPE GVCFS
// ------------------------------------------------------------------
process GENOTYPE_GVCFS {
    tag "$sample_id"
    publishDir "${params.outdir}/vcfs/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(gvcf), path(idx)
    path ref
    path ref_dict
    path ref_idx

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: raw_vcf

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \
        -R ${ref} \
        -V ${gvcf} \
        -O ${sample_id}.vcf.gz \
        -A AlleleFraction -A FisherStrand -A StrandOddsRatio -A ExcessHet
    """
}

// ------------------------------------------------------------------
// PROCESS 6: SPLIT AND FILTER
// ------------------------------------------------------------------
process FILTER_VARIANTS {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(idx)
    path ref
    path ref_dict
    path ref_idx

    output:
    tuple val(sample_id), path("${sample_id}.snps.filtered.vcf.gz"), path("${sample_id}.indels.filtered.vcf.gz"), emit: filtered_vcfs

    script:
    """
    # 1. Select SNPs
    gatk SelectVariants -R $ref -V $vcf -select-type SNP -O snps.vcf.gz
    
    # 2. Filter SNPs (Matches SNV_filter.sh)
    gatk VariantFiltration \
        -R $ref \
        -V snps.vcf.gz \
        -O ${sample_id}.snps.filtered.vcf.gz \
        --filter-name "QD2"             --filter-expression "QD < 2.0" \
        --filter-name "QUAL30"          --filter-expression "QUAL < 30.0" \
        --filter-name "SOR3"            --filter-expression "SOR > 3.0" \
        --filter-name "FS60"            --filter-expression "FS > 60.0" \
        --filter-name "MQ40"            --filter-expression "MQ < 40.0" \
        --filter-name "MQRankSum-12.5"  --filter-expression "MQRankSum < -12.5" \
        --filter-name "ReadPosRankSum-8" --filter-expression "ReadPosRankSum < -8.0"

    # 3. Select Indels
    gatk SelectVariants -R $ref -V $vcf -select-type INDEL -O indels.vcf.gz

    # 4. Filter Indels (Matches commented section in SNV_filter.sh)
    gatk VariantFiltration \
        -R $ref \
        -V indels.vcf.gz \
        -O ${sample_id}.indels.filtered.vcf.gz \
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

    FASTP(read_pairs_ch)
    ALIGN_AND_PREPROCESS(FASTP.out.trimmed_reads, ref_ch, ref_idx_ch, ref_dict_ch)
    HAPLOTYPE_CALLER(ALIGN_AND_PREPROCESS.out.bam_tuple, intervals_ch, ref_ch, ref_dict_ch, ref_idx_ch)
    MERGE_GVCFS(HAPLOTYPE_CALLER.out.gvcfs.groupTuple(), ref_dict_ch)
    GENOTYPE_GVCFS(MERGE_GVCFS.out.merged_gvcf, ref_ch, ref_dict_ch, ref_idx_ch)
    FILTER_VARIANTS(GENOTYPE_GVCFS.out.raw_vcf, ref_ch, ref_dict_ch, ref_idx_ch)
}
