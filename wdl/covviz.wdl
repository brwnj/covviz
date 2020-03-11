task run_indexcov {
    # input {
        File fasta_index
        Array[File] alignment_indexes
        String sexchroms = "X,Y"
        String excludepatt = "^GL|^hs|^chrEBV$|M$|MT$|^NC|_random$|Un_|^HLA\\-|_alt$|hap\\d+$"
        String project = "covviz"

        Int disk_size = 20
        Int memory = 8
    # }
    command {
        goleft indexcov --sex '${sexchroms}' ${"--excludepatt='" + excludepatt + "'"} \
            --directory ${project} --fai ${fasta_index} ${sep=" " alignment_indexes}
    }
    runtime {
        memory: memory + "GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        docker: "brwnj/covviz:v1.2.1"
    }
    output {
        Array[File] indexcov_pngs = glob("${project}/${project}-indexcov-*.png")
        Array[File] indexcov_html = glob("${project}/*.html")
        File indexcov_bed = "${project}/${project}-indexcov.bed.gz"
        File indexcov_ped = "${project}/${project}-indexcov.ped"
        File indexcov_roc = "${project}/${project}-indexcov.roc"
    }
    meta {
        author: "Joe Brown"
        email: "brwnjm@gmail.com"
        description: "Run @brentp's indexcov across alignment index"
    }
}

task run_covviz {
    # input {
        File bed
        File ped
        File? gff
        Int minsamples = 8
        String sexchroms = "X,Y"
        String excludepatt = "^GL|^hs|^chrEBV$|M$|MT$|^NC|_random$|Un_|^HLA\\-|_alt$|hap\\d+$"
        Float zthreshold = 3.5
        Int distancethreshold = 150000
        Int slop = 500000
        Boolean skipnorm = true

        Int disk_size = 20
        Int memory = 8
    # }
    command {
        covviz --min-samples ${minsamples} --sex-chroms '${sexchroms}' \
            ${"--exclude '" + excludepatt + "'"} --z-threshold ${zthreshold} \
            --distance-threshold ${distancethreshold} --slop ${slop} \
            --ped ${ped} ${"--gff " + gff} \
            ${true="--skip-norm" false="" skipnorm} ${bed}
    }
    runtime {
        memory: memory + "GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        docker: "brwnj/covviz:v1.2.1"
    }
    output {
        File covviz_report = "covviz_report.html"
    }
    parameter_meta {
        gff: "file path to gff matching genome build of `bed`"
        sexchroms: "sex chromosomes as they are in `bed`"
        excludepatt: "regular expression of chromosomes to skip"
        zthreshold: "a sample must greater than this many standard deviations in order to be found significant"
        distancethreshold: "consecutive significant points must span this distance in order to pass this filter"
        slop: "leading and trailing segments added to significant regions to make them more visible"
    }
    meta {
        author: "Joe Brown"
        email: "brwnjm@gmail.com"
        description: "Generate covviz HTML report from normalized input"
    }
}

workflow covviz {
    # input {
        # single column file with cram indexes in column 0
        File manifest
        Array[Array[String]] sample_data = read_tsv(manifest)
        File fasta_index

        String project = "covviz"
        File? gff
        Int minsamples = 8
        String sexchroms = "X,Y"
        String excludepatt = "^GL|^hs|^chrEBV$|M$|MT$|^NC|_random$|Un_|^HLA\\-|_alt$|hap\\d+$"
        Float zthreshold = 3.5
        Int distancethreshold = 150000
        Int slop = 500000

        Int disk_size = 20
        Int memory = 8
    # }
    call run_indexcov {
        input:
            fasta_index = fasta_index,
            alignment_indexes = flatten(sample_data),
            sexchroms = sexchroms,
            excludepatt = excludepatt,
            project = project,

            disk_size = disk_size,
            memory = memory
    }
    call run_covviz {
        input:
            bed = run_indexcov.indexcov_bed,
            ped = run_indexcov.indexcov_ped,
            gff = gff,
            minsamples = minsamples,
            sexchroms = sexchroms,
            excludepatt = excludepatt,
            zthreshold = zthreshold,
            distancethreshold = distancethreshold,
            slop = slop,

            disk_size = disk_size,
            memory = memory
    }
    meta {
        author: "Joe Brown"
        email: "brwnjm@gmail.com"
        description: "Multi-sample genome coverage viewer to observe large, coverage-based anomalies"
    }
}
