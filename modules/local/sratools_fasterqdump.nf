
process SRATOOLS_FASTERQDUMP {
    tag "$meta.id"
    label 'process_medium'

    container "ncbi/sra-tools:3.0.0"

    input:
    tuple val(meta), val(id)

    output:
    tuple val(meta), path(fastq_output), emit: reads
    tuple val(meta), path(md5_output)  , emit: md5
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"
    // Paired-end data extracted by fasterq-dump (--split-3 the default) always creates
    // *_1.fastq *_2.fastq files but sometimes also an additional *.fastq file
    // for unpaired reads which we ignore here.
    fastq_output = meta.single_end ? '*.fastq.gz'     : '*_{1,2}.fastq.gz'
    md5_output = meta.single_end ? '*.fastq.gz.md5' : '*_{1,2}.fastq.gz.md5'
    """
    eval "\$(vdb-config -o n NCBI_SETTINGS | sed 's/[" ]//g')"
    if [[ ! -f "\${NCBI_SETTINGS}" ]]; then
        mkdir -p "\$(dirname "\${NCBI_SETTINGS}")"
        printf '${config}' > "\${NCBI_SETTINGS}"
    fi

    fastq-dump \\
        $args \\
        ${id} \\
        --split-3 \\
        --unaligned

    pigz \\
        $args2 \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq

    ## Rename FastQ files by meta.id
    if [ -f  ${id}.fastq.gz ]; then
        mv ${id}.fastq.gz ${meta.id}.fastq.gz
        md5sum ${meta.id}.fastq.gz > ${meta.id}.fastq.gz.md5
    fi

    if [ -f  ${id}_1.fastq.gz ]; then
        mv ${id}_1.fastq.gz ${meta.id}_1.fastq.gz
        md5sum ${meta.id}_1.fastq.gz > ${meta.id}_1.fastq.gz.md5
    fi

    if [ -f  ${id}_2.fastq.gz ]; then
        mv ${id}_2.fastq.gz ${meta.id}_2.fastq.gz
        md5sum ${meta.id}_2.fastq.gz > ${meta.id}_2.fastq.gz.md5
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
