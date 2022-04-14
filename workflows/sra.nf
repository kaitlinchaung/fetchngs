/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def valid_params = [
    ena_metadata_fields : ['run_accession', 'experiment_accession', 'library_layout', 'fastq_ftp', 'fastq_md5']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSra.initialise(params, log, valid_params)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { SRP_TO_SRR              } from '../modules/local/srp_to_srr'
include { SRA_IDS_TO_RUNINFO      } from '../modules/local/sra_ids_to_runinfo'
include { SRA_RUNINFO_TO_FTP      } from '../modules/local/sra_runinfo_to_ftp'
include { SRA_FASTQ_FTP           } from '../modules/local/sra_fastq_ftp'
include { SRA_TO_SAMPLESHEET      } from '../modules/local/sra_to_samplesheet'
include { SRA_MERGE_SAMPLESHEET   } from '../modules/local/sra_merge_samplesheet'
include { MULTIQC_MAPPINGS_CONFIG } from '../modules/local/multiqc_mappings_config'

include { SRA_FASTQ_SRATOOLS      } from '../subworkflows/local/sra_fastq_sratools'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SRA {


    main:

    // Read in inputs to ch_ids
    if (params.srp) {
        //
        // MODULE: Get SRR numbers from SRP project
        //
        SRP_TO_SRR (
            params.srp
        )

        id_list = SRP_TO_SRR.out.ids

    } else {
        // use input id list
        id_list = Channel.fromPath(params.id_list)

    }

    // Read in ids from SRR numbers
    id_list
        .splitCsv(
            header: false,
            sep: '',
            strip:true
        )
        .map { it[0] }
        .set { ids }

    if (params.num_samples) {
        ids = ids.take(params.num_samples)
    }

    //
    // MODULE: Get SRA run information for public database ids
    //
    SRA_IDS_TO_RUNINFO (
        ids,
        params.ena_metadata_fields ?: ''
    )

    //
    // MODULE: Parse SRA run information, create file containing FTP links and read into workflow as [ meta, [reads] ]
    //
    SRA_RUNINFO_TO_FTP (
        SRA_IDS_TO_RUNINFO.out.tsv
    )

    // Concatenate all metadata files into 1 mega file
    SRA_RUNINFO_TO_FTP.out.tsv
        .map { file ->
            file.text + '\n'
        }
        .collectFile (
            name:       "metadata.tsv",
            storeDir:   "${params.outdir}/metadata",
            keepHeader: true,
            skip:       1
        )

    SRA_RUNINFO_TO_FTP
        .out
        .tsv
        .splitCsv(header:true, sep:'\t')
        .map {
            meta ->
                meta.single_end = meta.single_end.toBoolean()
                [ meta, [ meta.fastq_1, meta.fastq_2 ] ]
        }
        .unique()
        .branch {
            ftp: it[0].fastq_1  && !params.force_sratools_download
            sra: !it[0].fastq_1 || params.force_sratools_download
        }
        .set { ch_sra_reads }


    if (!params.skip_fastq_download) {

        if (params.unaligned) {
            //
            // SUBWORKFLOW: Download sequencing reads without FTP links using sra-tools.
            //
            SRA_FASTQ_SRATOOLS (
                ch_sra_reads.sra.map { meta, reads -> [ meta, meta.run_accession ] }
            )

        } else {
            //
            // MODULE: If FTP link is provided in run information then download FastQ directly via FTP and validate with md5sums
            //
            SRA_FASTQ_FTP (
                ch_sra_reads.ftp
            )
        }

    }

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    WorkflowSra.curateSamplesheetWarn(log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
