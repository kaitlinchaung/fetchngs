
process SRP_TO_SRR {
    label 'error_retry'

    conda (params.enable_conda ? "conda-forge::python=3.9.7 bioconda::pysradb" : null)

    input:
    val srp

    output:
    path id_list, emit: ids

    script:
    id_list = "ids_${srp}.txt"
    """
    pysradb srp-to-srr ${srp} | awk '(NR>1) {print \$2}' > ${id_list}
    """
}
