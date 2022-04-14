//
// Download FASTQ sequencing reads from the NCBI's Sequence Read Archive (SRA).
//

include { SRATOOLS_PREFETCH    } from '../../modules/local/sratools_prefetch'
include { SRATOOLS_FASTERQDUMP } from '../../modules/local/sratools_fasterqdump'

workflow SRA_FASTQ_SRATOOLS {
    take:
    sra_ids  // channel: [ val(meta), val(id) ]

    main:

    //
    // Prefetch sequencing reads in SRA format.
    //
    SRATOOLS_PREFETCH ( sra_ids )

    //
    // Convert the SRA format into one or more compressed FASTQ files.
    //
    SRATOOLS_FASTERQDUMP ( SRATOOLS_PREFETCH.out.sra )

    emit:
    reads    = SRATOOLS_FASTERQDUMP.out.reads  // channel: [ val(meta), [ reads ] ]
}
