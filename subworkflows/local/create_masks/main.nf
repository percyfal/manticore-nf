include { BEDTOOLS_MASKFASTA as BEDTOOLS_MASKFASTA_UNMASKED_TEMPLATE } from '../../../modules/local/bedtools/maskfasta/main'
include { BEDTOOLS_MASKFASTA as BEDTOOLS_MASKFASTA_MASKED_TEMPLATE } from '../../../modules/local/bedtools/maskfasta/main'
include { BEDTOOLS_MASKFASTA as BEDTOOLS_MASKFASTA_UNMASK } from '../../../modules/local/bedtools/maskfasta/main'
include { BEDTOOLS_MASKFASTA as BEDTOOLS_MASKFASTA_MASK } from '../../../modules/local/bedtools/maskfasta/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_CONVERT_MULTILINE_FASTA } from '../../../modules/local/seqkit/seq/main'

workflow CREATE_MASKS {

    take:
    genome_bed         // channel: meta, bed
    fasta              // channel: fasta
    mask_bed           // channel: meta, bed
    windowed_mask_bed           // channel: meta, bed - FIXME: is this
				// even needed?

    main:

    ch_versions = Channel.empty()

    // Make template mask file with all positions masked
    BEDTOOLS_MASKFASTA_MASKED_TEMPLATE(genome_bed, fasta)
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA_MASKED_TEMPLATE.out.versions)

    // Generate mask for both coverage and without
    BEDTOOLS_MASKFASTA_UNMASK(mask_bed, BEDTOOLS_MASKFASTA_MASKED_TEMPLATE.out.fasta.map{it[1]})

    // Somehow collect all different variants of mask files in two
    // channels: one for windowed analyses, one for non-windowed
    // analyses; these are emitted below

    // convert multiline fasta to one-line
    SEQKIT_SEQ_CONVERT_MULTILINE_FASTA(BEDTOOLS_MASKFASTA_UNMASK.out.fasta)

    // Partition bed file by feature class


    emit:
    fasta      = SEQKIT_SEQ_CONVERT_MULTILINE_FASTA.out.fasta           // channel: [ val(meta), [ fasta ] ]
    windowed_fasta = SEQKIT_SEQ_CONVERT_MULTILINE_FASTA.out.fasta           // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions                                     // channel: [ versions.yml ]
}
