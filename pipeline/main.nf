/*  CNV Detection Pipeline  
 *  Usage:  
 *
 *  Author: Gunjun Sethia
 *  Modified by Angela Hickey
 *  NYU Center for Genetics and System Biology 2020
 */

// set some parameters
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/cnv_temp"
params.minimum_read_length = 50

// Define modules here
BWA = 'bwa/intel/0.7.15'
TRIMGALORE = 'trim_galore/0.4.4'
CUTADAPT = 'cutadapt/intel/1.12'
SAMTOOLS = 'samtools/intel/1.3.1'
PICARD_JAR = '/share/apps/picard/2.8.2/picard-2.8.2.jar'


// print starting parameters
println "reads: $params.reads"
println "ref: $params.ref"
println "output: $params.out"
println "cnv temp dir: $params.tmpdir"

// Setup the reference file
reference_ch = Channel.fromPath(params.ref)
reference_ch.into { ref_align; ref_insert_metric }

/* Prepare the fastq read pairs for input.
 * While doing this, count number of input samples
 */
num_samples = 0
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .tap { read_pairs_ch }
    .subscribe({ num_samples += 1 })

process trim {
    publishDir "${params.out}/trimmed_reads", mode:'copy'

    input:
    set pair_id, file(reads) from read_pairs_ch

    output:
    set val(pair_id), file("${pair_id}01_AP1-A02_val_1.fq.gz"), file("${pair_id}02_AP1-A02_val_2.fq.gz") into trimmed_reads_ch

    script:
    """
    module purge
    module load $TRIMGALORE
    module load $CUTADAPT
    module load fastqc/0.11.5

    trim_galore --fastqc --length ${params.minimum_read_length} \
    --trim-n \
    --clip_R1 ${params.bases_trim_5prime_Read1} \
    --clip_R2 ${params.bases_trim_5prime_Read2} \
    --three_prime_clip_R1 ${params.bases_trim_3prime_Read1} \
    --three_prime_clip_R2 ${params.bases_trim_3prime_Read2} \
    --paired \
    ${reads[0]} \
    ${reads[1]}

    """
}    

process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'
	
    input:
    set pair_id, file(trimmed_reads1), file(trimmed_reads2) from trimmed_reads_ch 
    file (ref) from ref_align
     
    output:
    set val(pair_id), file("${pair_id}_aligned_reads.sorted.bam"), file("${pair_id}_aligned_reads.sorted.bam.bai") into aligned_reads_ch

    script:
    """
    module purge
    module load $BWA
    module load $SAMTOOLS
    bwa index ${ref}
    bwa mem -t 20 ${ref} ${trimmed_reads1} ${trimmed_reads2} > ${pair_id}_aligned_reads.bam
    samtools sort ${pair_id}_aligned_reads.bam > ${pair_id}_aligned_reads.sorted.bam
    samtools index ${pair_id}_aligned_reads.sorted.bam
    
    """
}

process insert_size {
    publishDir "${params.out}/insert_size_metrics", mode:'copy'

    input:
    set pair_id, file(sorted_bam), file(index_bam) from aligned_reads_ch
    file (ref) from ref_insert_metric 

    output:
    set val(pair_id), file("${pair_id}_alignment_metrics.txt"), file("${pair_id}_insert_size_histogram.pdf"), file("${pair_id}_insert_size_metrics.txt") into insert_size_ch

    script: 
    """
    # obtaining alignment metrics using Picards tools
    module purge
    module load picard/2.8.2
    module load r/intel/3.3.2
    java -jar $PICARD_JAR \
    CollectAlignmentSummaryMetrics \
    R=$ref \
    I=$sorted_bam \
    O=alignment_metrics.txt

    grep -v '^#' alignment_metrics.txt | cut -f1-12 | sed '/^[[:space:]]*\$/d' > ${pair_id}_alignment_metrics.txt

    # obtaining insert size metrics using Picards tools
    java -jar $PICARD_JAR \
    CollectInsertSizeMetrics \
    INPUT=$sorted_bam \
    OUTPUT=insert_metrics.txt \
    HISTOGRAM_FILE=${pair_id}_insert_size_histogram.pdf

    head -n 9 insert_metrics.txt | grep -v '^#' | cut -f1-19|sed '/^[[:space:]]*\$/d' >${pair_id}_insert_size_metrics.txt
    
    """
}
