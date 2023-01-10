include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                         } from '../../../modules/nf-core/samtools/faidx/main'
include { BUILD_INTERVALS as BUILD_INTERVALS_FAIDX      } from '../../../modules/local/build_intervals/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_GENOME_BED } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow PREPARE_GENOME {
    take:
        fasta                   // channel: [mandatory] fasta
        fasta_fai               // channel: [optional]  fasta_fai

    main:

    ch_versions = Channel.empty()

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].baseName], it] })

    fasta_fai = SAMTOOLS_FAIDX.out.fai
    BUILD_INTERVALS_FAIDX(fasta_fai.map{ it -> [[id: "${it[0].id}.genome"], it[1]]})

    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(BUILD_INTERVALS_FAIDX.out.versions)

    emit:
    dict                             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
    fasta_fai                        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                      // path: genome.fasta.fai
    genome_bed                       = BUILD_INTERVALS_FAIDX.out.bed                        // path: genome.fasta.bed
    genome_txt                       = BUILD_INTERVALS_FAIDX.out.txt    // path: genome.fasta.txt
    versions                         = ch_versions
}
