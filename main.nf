#!/usr/bin/env nextflow

// required arguments
params.indexes = false
if( !params.indexes ) { exit 1, "--indexes is not defined" }
params.fai = false
if( !params.fai ) { exit 1, "--fai is not defined" }
params.gff = false
if (!params.gff ) { exit 1, "--gff is not defined" }
params.ped = false

project = params.project ?: 'NF'

log.info("\n")
log.info("====================================================================")
log.info("covviz - find large, coverage-based variations on chromosomes ")
log.info("====================================================================")
log.info("\n")
log.info("#### Homepage / Documentation")
log.info("https://github.com/brwnj/covviz")
log.info("#### Authors")
log.info("Joe Brown <brwnjm@gmail.com>")
log.info("\n")
log.info("====================================================================")
log.info("\n")
log.info("Alignment indexes            --indexes          : ${params.indexes}")
log.info("Reference index              --fai              : ${params.fai}")
log.info("GFF                          --gff              : ${params.gff}")
log.info("Sex chromosomes              --sexchroms        : ${params.sexchroms}")
log.info("Excluded chroms              --exclude          : ${params.exclude}")
log.info("Z threshold                  --zthreshold       : ${params.zthreshold}")
log.info("Distance threshold           --distancethreshold: ${params.distancethreshold}")
log.info("Significant region slop      --slop             : ${params.slop}")
log.info("\n")
log.info("====================================================================")
log.info("covviz is free and unrestricted for non-commercial use.       ")
log.info("For commercial use, please contact [bpedersen@base2genomics.com].  ")
log.info("====================================================================")
log.info("\n")

// instantiate files
fai = file(params.fai)
gff = file(params.gff)
outdir = file(params.outdir)
custom_ped = false
if (params.ped) {
    custom_ped = file(params.ped)
}

// check file existence
if( !fai.exists() ) { exit 1, "Reference FASTA [${fai}] index does not exist." }

Channel
    .fromPath(params.indexes, checkIfExists: true)
    .set { index_ch }

process run_indexcov {
    publishDir path: "$outdir/indexcov", mode: "copy"
    container 'brentp/smoove:v0.2.3'
    memory 8.GB
    cache 'deep'

    input:
    file idx from index_ch.collect()
    file fai

    output:
    file("${project}*.png")
    file("*.html")
    file("${project}*.bed.gz") into bed_ch
    file("${project}*.ped") into ped_ch
    file("${project}*.roc") into roc_ch

    script:
    excludepatt = params.exclude ? "--excludepatt \"${params.exclude}\"" : ''
    """
    goleft indexcov --sex ${params.sexchroms} $excludepatt --directory $project --fai $fai $idx
    mv $project/* .
    """
}

// account for optional, custom ped and the need to merge that with indexcov output
(merge_ch, report_ch) = (params.ped ? [ped_ch, Channel.empty()]: [Channel.empty(), ped_ch])

process merge_peds {
    input:
    file ped from merge_ch
    file custom_ped

    output:
    file 'merged.ped' into merged_ch

    script:
    template 'merge_peds.py'
}

process build_report {
    publishDir path: "$outdir", mode: "copy", pattern: "*.html", overwrite: true
    container 'brwnj/covviz:v1.0.5'

    input:
    file ped from report_ch.mix(merged_ch).collect()
    file roc from roc_ch
    file bed from bed_ch
    file gff
    file custom_ped

    output:
    file("covviz_report.html")

    script:
    """
    covviz --min-samples ${params.minsamples} --sex-chroms ${params.sexchroms} --exclude '${params.exclude}' \
        --z-threshold ${params.zthreshold} --distance-threshold ${params.distancethreshold} \
        --slop ${params.slop} --ped ${ped} --gff ${gff} ${bed} ${roc}
    """
}
