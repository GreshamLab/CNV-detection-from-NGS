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
    module purge
    module load trim_galore/0.4.4
    module load cutadapt/intel/1.12

    trim_galore --fastqc --length ${minimum_read_length} \
    --trim-n \
    --clip_R1 ${bases_trim_5prime_Read1} \
    --clip_R2 ${bases_trim_5prime_Read2} \
    --three_prime_clip_R1 ${bases_trim_3prime_Read1} \
    --three_prime_clip_R2 ${bases_trim_3prime_Read2} \
    -o $PWD \
    --paired \
    ${fastq1} \
    ${fastq2}

}
