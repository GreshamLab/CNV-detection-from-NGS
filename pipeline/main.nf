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
reference_ch.into { ref_align; ref_insert_metric; ref_pindel }

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

// Copy aligned_reads_ch to be used in multiple processes
( bam_insert_size, bam_read_depth, bam_pindel, bam_lumpy ) = aligned_reads_ch.into(4)

process insert_size {
    publishDir "${params.out}/insert_size_metrics", mode:'copy'

    input:
    set pair_id, file(sorted_bam), file(index_bam) from bam_insert_size
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

process read_depth {
    publishDir "${params.out}/read_depth", mode:'copy'

    input:
    set pair_id, file(sorted_bam), file(index_bam) from bam_read_depth

    output: 
    set val(pair_id), file("${pair_id}_RD.txt") into read_depth_ch

    script:
    """
    # obtaining read depth ie coverage using samtools
    module load samtools/intel/1.3.1
    samtools depth -a $sorted_bam > ${pair_id}_RD.txt
    """
}

process pindel {
    publishDir "${params.out}/pindel", mode:'copy'

    input:
    set pair_id, file(sorted_bam), file(index_bam) from bam_pindel
    file(ref) from ref_pindel
    file(alignment_metrics), file(histogram), file(insert_size_metrics) from (insert_size_ch)

    output:
    set val(pair_id), file("${pair_id}_pindel.vcf") into pindel_ch

    script:
    """
    # obtain config file for pindel
    mean_IS = /$(sed -n '2p' < $insert_size_metrics | cut -f 5)
    echo "$sorted_bam ${mean_IS} ${pair_id}" > config_${pair_id}.txt

    module purge
    module load pindel/intel/20170402
    /share/apps/pindel/20170402/intel/bin/pindel \
    -T 20 \
    -f $ref \
    -i config_${pair_id}.txt \
    -c ALL \
    -o ${pair_id}_output

    pindel2vcf -p ${pair_id}_output_D -r $ref -R UCSC_SacCer -d Feb2017 -v ${pair_id}_DEL_pindel.vcf
    pindel2vcf -p ${pair_id}_output_TD -r $ref -R UCSC_SacCer -d Feb2017 -v ${pair_id}_DUP_pindel.vcf

    module purge
    module load vcftools/intel/0.1.14
    vcf-concat ${pair_id}_DEL_pindel.vcf ${pair_id}_DUP_pindel.vcf > ${pair_id}_pindel.vcf

    echo "pindel done"
    
    """
}