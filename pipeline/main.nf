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


// print starting parameters
println "reads: $params.reads"
println "ref: $params.ref"
println "output: $params.out"
println "cnv temp dir: $params.tmpdir"

// Setup the reference file
ref = file(params.ref)

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
    set val(pair_id), file("${pair_id}_val_1.fq.gz"), file("${pair_id}_val_2.fq.gz) \
	into trimmed_reads_ch

    script:
    """
    module purge
    module load $TRIMGALORE
    module load $CUTADAPT

    trim_galore --fastqc --length ${params.minimum_read_length} \
    --trim-n \
    --clip_R1 ${params.bases_trim_5prime_Read1} \
    --clip_R2 ${params.bases_trim_5prime_Read2} \
    --three_prime_clip_R1 ${params.bases_trim_3prime_Read1} \
    --three_prime_clip_R2 ${params.bases_trim_3prime_Read2} \
    -o $PWD \
    --paired \
    ${reads[0]} \        
    ${reads[1]}

    """

process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'
	
    input:
    set pair_id, file(trimmed_reads) from trimmed_reads_ch
     
    output:
    set val(pair_id), file("${pair_id}_aligned_reads.bam") \
	into aligned_reads_ch

    script:
    """
    module purge
    module load $BWA
    module load $SAMTOOLS
    bwa mem -t 20 ${ref} ${read[0]} ${read[1]} > ${pair_id}.bam
    samtools sort ${pair_id}.bam > ${pair_id}.sorted.bam
    samtools index ${pair_id}.sorted.bam
    
    """


}
