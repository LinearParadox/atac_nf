include { build_index } from '../modules/bowtie2/align.nf'
include { align } from '../modules/bowtie2/align.nf'

workflow qc_samples {
    take:
        samples
        index
        fasta
    main:
    if(file(index).exists()){
        if (index.isDirectory){
            index = channel.fromPath(index+"/*bt2").collect()
        }
    } else{
        index = build_index(fasta).collect()
    }
    aligned = align(samples, index)
    emit:

    }