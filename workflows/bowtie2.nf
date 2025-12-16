include { build_index } from '../modules/bowtie2/align.nf'
include { align } from '../modules/bowtie2/align.nf'

workflow bowtie2 {
    take:
        samples
        index
        fasta
    main:
    if(file(index).exists()){
        if (index.isDirectory){
            index = channel.fromPath(index+"/*bt2").collect()
        } else{
            index = channel.fromPath(index).collect()
        }
    } else{
        index = build_index(fasta).collect()
    }
    aligned = align(samples, params.fragment_size, params.multimap, index)


    emit:
        aligned_bam = aligned.aligned_bam
        alignment_metrics = aligned.alignment_metrics
    }