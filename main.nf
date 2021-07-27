#!/usr/bin/env nextflow

params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------

    Required arguments

    --indexes              quoted file path with wildcard ('*.crai') to
                           cram or bam indexes
    --fai                  file path to .fai reference index


    Workflow Options

    --outdir               output directory for results
                           default: "./results"
    --gff                  file path to gff matching genome build of
                           `--indexes`
    --sexchroms            sex chromosomes as they are in `--indexes`
                           default: "X,Y"
    --exclude              regular expression of chromosomes to skip
                           default: "^GL|^hs|^chrEBV\$|M\$|MT\$|^NC|_random\$|Un_|^HLA\\-|_alt\$|hap\\d+\$"
    --zthreshold           a sample must greater than this many standard
                           deviations in order to be found significant
                           default: 3.5
    --distancethreshold    consecutive significant points must span this
                           distance in order to pass this filter
                           default: 150000
    --slop                 leading and trailing segments added to
                           significant regions to make them more visible
                           default: 500000
    --ped                  custom metadata that will be merged with the
                           .ped output of indexcov
                           default: false
    --samplecol            the column header for sample IDs in your custom
                           ped file
                           default: "sample_id"

    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// required arguments
params.indexes = false
if( !params.indexes ) { exit 1, "--indexes is not defined" }
params.fai = false
if( !params.fai ) { exit 1, "--fai is not defined" }
params.gff = false
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
if (params.gff) {
log.info("GFF                          --gff              : ${params.gff}")
}
log.info("Sex chromosomes              --sexchroms        : ${params.sexchroms}")
log.info("Excluded chroms              --exclude          : ${params.exclude}")
log.info("Z threshold                  --zthreshold       : ${params.zthreshold}")
log.info("Distance threshold           --distancethreshold: ${params.distancethreshold}")
log.info("Significant region slop      --slop             : ${params.slop}")
log.info("\n")
log.info("====================================================================")
log.info("\n")

// instantiate files
fai = file(params.fai)
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

Channel
    .from(params.gff ? file(params.gff) : false)
    .set { gff_ch }

process run_indexcov {
    publishDir path: "$outdir/indexcov", mode: "copy"
    label 'indexcov'

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

    input:
    file ped from report_ch.mix(merged_ch).collect()
    file roc from roc_ch
    file bed from bed_ch
    file gff from gff_ch

    output:
    file("covviz_report.html")

    script:
    gff_opt = params.gff ? "--gff ${gff}" : ""
    gff_attr = params.gffattr ? "--gff-attr ${gffattr}" : ""
    """
    covviz --min-samples ${params.minsamples} --sex-chroms ${params.sexchroms} --exclude '${params.exclude}' \
        --z-threshold ${params.zthreshold} --distance-threshold ${params.distancethreshold} \
        --slop ${params.slop} --ped ${ped} ${gff_opt} ${gff_attr} --skip-norm ${bed}
    """
}
