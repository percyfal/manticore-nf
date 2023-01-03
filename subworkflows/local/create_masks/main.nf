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

    main:

    ch_versions = Channel.empty()

    // Add genome_bed to mask_bed list
    mask_bed_combined = mask_bed.mix(genome_bed.map{[[id: "genome", mode: "site"], it[1]]}).mix(genome_bed.map{[[id: "genome", mode: "window"], it[1]]})
    mask_bed_grouped = mask_bed_combined.groupTuple(by: 1).map{[[id: it[0][0].id, mode: it[0].mode, bed: it[1]], it[1]]}

    // Make template mask file with all positions masked
    BEDTOOLS_MASKFASTA_MASKED_TEMPLATE(genome_bed, fasta)
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA_MASKED_TEMPLATE.out.versions)

    BEDTOOLS_MASKFASTA_UNMASK(mask_bed_grouped, BEDTOOLS_MASKFASTA_MASKED_TEMPLATE.out.fasta.map{it[1]})

    // convert multiline fasta to one-line
    SEQKIT_SEQ_CONVERT_MULTILINE_FASTA(BEDTOOLS_MASKFASTA_UNMASK.out.fasta)
    ch_versions = ch_versions.mix(SEQKIT_SEQ_CONVERT_MULTILINE_FASTA.out.versions)


    emit:
    fasta    = SEQKIT_SEQ_CONVERT_MULTILINE_FASTA.out.fasta           // channel: [ val(meta), fasta  ]
    versions = ch_versions                                              // channel: [ versions.yml ]
}
