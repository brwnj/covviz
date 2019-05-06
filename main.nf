#!/usr/bin/env nextflow

// required arguments
params.indexes = false
if( !params.indexes ) { exit 1, "--indexes is not defined" }
params.fai = false
if( !params.fai ) { exit 1, "--fai is not defined" }
params.gff = false
if (!params.gff ) { exit 1, "--gff is not defined" }

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

// check file existence
if( !fai.exists() ) { exit 1, "Reference FASTA [${fai}] index does not exist." }

Channel
    .fromPath(params.indexes, checkIfExists: true)
    .set { index_ch }

process run_indexcov {
    publishDir path: "$outdir/indexcov", mode: "copy"
    container 'brentp/smoove:v0.2.3'
    memory 4.GB
    cache 'deep'

    input:
    file idx from index_ch.collect()
    file fai

    output:
    file("${project}*.png")
    file("*.html")
    file("${project}*.bed.gz") into coverage_bed
    file("${project}*.ped") into ped
    file("${project}*.roc") into roc

    script:
    excludepatt = params.exclude ? "--excludepatt \"${params.exclude}\"" : ''
    """
    goleft indexcov --sex ${params.sexchroms} $excludepatt --directory $project --fai $fai $idx
    mv $project/* .
    """
}

process build_report {
    publishDir path: "$outdir", mode: "copy", pattern: "*.html", overwrite: true
    container 'jupyter/scipy-notebook:7d427e7a4dde'

    input:
    file pedfile from ped
    file rocfile from roc
    file bedfile from coverage_bed
    file gff

    output:
    file("covviz_report.html")

    script:
    template 'parse_indexcov.py'
}
