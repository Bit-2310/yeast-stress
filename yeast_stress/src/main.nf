nextflow.enable.dsl=2

params.data_dir = './data'
params.results_dir = './results'

process FastQC {
    container 'fastqc'
    input:
    path reads

    output:
    path "${reads.baseName}_fastqc.html", emit: report
    path "${reads.baseName}_fastqc.zip", emit: qc

    script:
    """
    fastqc ${reads} -o ./
    """
}

process HISAT2_Align {
    container 'hisat2'
    input:
    tuple val(id), path(reads)
    path genome_index

    output:
    path "${id}_aligned.sam", emit: aligned

    script:
    """
    hisat2 -x ${genome_index} -U ${reads} -S ${id}_aligned.sam
    """
}

process SortBAM {
    container 'samtools'
    input:
    path sam_file

    output:
    path "${sam_file.baseName}_sorted.bam", emit: sorted

    script:
    """
    samtools sort ${sam_file} -o ${sam_file.baseName}_sorted.bam
    """
}

workflow {
    data = file("${params.data_dir}/*.fastq")
    genome_index = file("${params.data_dir}/genome/sacCer3")

    fastqc_results = FastQC(data)
    aligned_sam = HISAT2_Align(data, genome_index)
    sorted_bam = SortBAM(aligned_sam)
}
