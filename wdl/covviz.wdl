task run_indexcov {
    input {
        File fasta_index
        Array[File] alignment_indexes
        String sexchroms = "X,Y"
        String excludepatt = "^GL|^hs|^chrEBV$|M$|MT$|^NC|_random$|Un_|^HLA\\-|_alt$|hap\\d+$"
        String project = "covviz"

        Int disk_size = 20
        Int memory = 8
    }
    command {
        goleft indexcov --sex ${sexchroms} ${true="--excludepatt" false="" excludepatt} ${excludepatt} \
            --directory ${project} --fai ${fasta_index} ${sep=" " alignment_indexes}
        mv ${project}/* .
    }
    runtime {
        memory: memory + "GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        docker: "brentp/smoove:v0.2.5"
    }
    output {
        Array[File] indexcov_pngs = glob("${project}-indexcov-*.png")
        Array[File] indexcov_html = glob("*.html")
        File indexcov_bed = "${project}-indexcov.bed.gz"
        File indexcov_ped = "${project}-indexcov.ped"
        File indexcov_roc = "${project}-indexcov.roc"
    }
    meta {
        author: "Joe Brown"
        email: "brwnjm@gmail.com"
        description: "Run @brentp's indexcov across alignment index"
    }
}

task run_covviz {
    input {
        File bed
        File ped
        File? gff
        Int minsamples = 8
        String sexchroms = "X,Y"
        String excludepatt = "^GL|^hs|^chrEBV$|M$|MT$|^NC|_random$|Un_|^HLA\\-|_alt$|hap\\d+$"
        Float zthreshold = 3.5
        Int distancethreshold = 150000
        Int slop = 500000
        Boolean skipnorm = "true"

        Int disk_size = 20
        Int memory = 8
    }
    command {
        covviz --min-samples ${minsamples} --sex-chroms ${sexchroms} \
            ${true="--excludepatt" false="" excludepatt} ${excludepatt} \
            --z-threshold ${zthreshold} --distance-threshold ${distancethreshold} \
            --slop ${slop} --ped ${ped} ${true="--gff" false="" gff} ${gff} \
            ${true="--skip-norm" false="" skipnorm} ${bed}
    }
    runtime {
        memory: memory + "GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        docker: "brwnj/covviz:v1.1.5"
    }
    output {
        File covviz_report = "covviz_report.html"
    }
    meta {
        author: "Joe Brown"
        email: "brwnjm@gmail.com"
        description: "Generate covviz HTML report from normalized"
    }
}

workflow smoove {

    File manifest
    Array[Array[File]] sample_data = read_tsv(manifest)
    File fasta
    File fasta_index
    File bed
    File gff
    File ? ped
    File ? known_sites
    String project_id
    String exclude_chroms
    Boolean sensitive

    Int small_disk
    Int medium_disk
    Int large_disk

    scatter (sample in sample_data) {
        call smoove_call {
            input:
                sample_id = sample[0],
                alignments = sample[1],
                alignments_index = sample[2],
                fasta = fasta,
                fasta_index = fasta_index,
                bed = bed,
                exclude_chroms = exclude_chroms,
                sensitive = sensitive,
                disk_size = small_disk
        }
    }
    call smoove_merge {
        input:
            project_id = project_id,
            fasta = fasta,
            fasta_index = fasta_index,
            vcfs = smoove_call.called_vcf,
            vcf_indexes = smoove_call.called_vcf_index,
            disk_size = small_disk
    }
    scatter (sample in sample_data) {
        call smoove_genotype {
            input:
                sample_id = sample[0],
                alignments = sample[1],
                alignments_index = sample[2],
                fasta = fasta,
                fasta_index = fasta_index,
                vcf = smoove_merge.merged_vcf,
                disk_size = small_disk
        }
    }
    call smoove_square {
        input:
            project_id = project_id,
            vcfs = smoove_genotype.genotyped_vcf,
            vcf_indexes = smoove_genotype.genotyped_vcf_index,
            gff = gff,
            disk_size = small_disk
    }
}


workflow {
    Array[File]
    File fasta_index

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

}
