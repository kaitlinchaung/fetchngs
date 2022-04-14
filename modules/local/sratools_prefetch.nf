
process SRATOOLS_PREFETCH {
    tag "$id"
    label 'process_low'
    label 'error_retry'

    container "ncbi/sra-tools:3.0.0"

    input:
    tuple val(meta), val(id)

    output:
    tuple val(meta), val(id)    , emit: sra
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"
    """
    eval "\$(vdb-config -o n NCBI_SETTINGS | sed 's/[" ]//g')"
    if [[ ! -f "\${NCBI_SETTINGS}" ]]; then
        mkdir -p "\$(dirname "\${NCBI_SETTINGS}")"
        printf '${config}' > "\${NCBI_SETTINGS}"
    fi

    prefetch \\
        --force all \\
        $args \\
        $id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}
