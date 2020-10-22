# covviz

Coverage visualization; a many-sample coverage browser.

The aim of `covviz` is to highlight regions of significant
(passing the user's z-score threshold) and sustained (beyond user specified
distance) deviation of coverage depth from the majority of samples. Significance is determined
using z-scores for all samples at all points using median absolute deviation.
In order for regions to be highlighted, points must be significant
consecutively throughout a user specified distance.

If you are analyzing a low number of samples, deviation may be irrelevant. In
this case, we can set `--min-samples` to be greater than our sample total
to skip Z-threshold calculation and plot coverages for all samples at all
points.

# Getting started

## From alignments (.bam and/or .cram)

Alignments must be indexed. The input for the `covviz` workflow are the indexes
of the alignments. For BAM, that would be .bai, and .crai for CRAM. Indexes
can be generated using [samtools](https://github.com/samtools/samtools) on your
sorted alignments:

```
samtools index mybam.bam
# generates mybam.bam.bai
```

### Installation and usage

Install [Nextflow](https://www.nextflow.io/) if you don't already have it. The only
dependency is Java 8 or later, then you run:

```
curl -s https://get.nextflow.io | bash
```

Or via [Bioconda](https://bioconda.github.io/recipes/nextflow/README.html) using:

```
conda install -c bioconda nextflow
```

Full nextflow installation instructions are available at:
https://www.nextflow.io/

There is no need to download the covviz code prior to execution or any software dependencies
when using a container service like Docker or Singularity.

### Docker/Singularity

To simplify prerequisite software installations and software version tracking,
we strongly recommend running `covviz` using Docker or Singularity. Docker
installation instructions for your operating system are available at:
https://docs.docker.com/install/

Then, with Docker or Singularity we run:

```
nextflow run brwnj/covviz -latest -profile docker \
    --indexes 'data/indexes/*.crai' \
    --fai data/g1k_v37_decoy.fa.fai \
    --gff data/Homo_sapiens.GRCh37.82.gff3.gz
```

Which gives us `./results/covviz_report.html`.

### Required arguments

+ `--indexes`
    + quoted file path with wildcard ('*.crai') to cram or bam indexes
+ `--fai`
    + file path to .fai reference index

A complete list of arguments can be displayed using:

```
nextflow run brwnj/covviz -latest --help
```

### Nextflow arguments

In the example above `-latest` pulls whatever the latest `covviz` code exists on GitHub
prior to execution and `-profile docker` sets `-with-docker` within Nextflow.

Other notable options are `-resume`, which when running a workflow a second will start
where previous runs of the workflow left off; and `-work-dir` which sets the location of
all intermediate files generated throughout the workflow.

## From coverage intervals (.bed)

The `covviz` CLI accepts bed3+ as input. If you've already generated your coverage
files you can start here and not the Nextflow workflow.

If you would prefer to run `indexcov` yourself across your .bai or .crai files,
the workflow above simply runs:

```
fai=data/g1k_v37_decoy.fa.fai
goleft indexcov --directory myproject --fai $fai *.crai
```

This will generate the expected inputs in their anticipated formats for the `covviz` CLI.

### Expected file format

To analyze your coverage data it needs to be in bed3+ format and include a
header with sample IDs. The first three column headers are agnostic, but
for samples test_sample1, test_sample2, and test_sample3, this would look like:

```
#chrom   start   end   sample1   sample2   sample3
```

### Installation of CLI and usage

To install the `covviz` Python package use:

```
pip install -U covviz
```

Then CLI usage is:

```
covviz $bed
```

A complete list of arguments can be displayed using:

```
covviz --help
```

### Adding custom metadata (.ped)

There is support for non-indexcov .ped files, though you may have to change
the default column IDs pertaining to the column which contains the sample ID
and the sex of the sample.

```
covviz --ped $ped --sample-col sample_col --sex sex_col $bed
```

### Adding annotation tracks

![significant_regions](data/img/covviz_tracks.gif)

Currently we support GFF, VCF, and BED. GFF tracks are added using `--gff`
where features are 'gene' and attributes have 'Name='. Feature type and
attribute regex can be configured using `--gff-feature` and `--gff-attr`.

VCF tracks (v4.1) are added with `--vcf` with the entire INFO string
being displayed by default. Specifying `--vcf-info` with something like
'CLNDN=' will grab just that field when using ClinVar variants. Including
large INFO strings for all variants can dramatically increase the size
of the covviz report.

Region based annotation tracks can be added using `--bed`. The name field
will be used to identify the regions when present.

Annotation tracks, `--gff`, `--vcf`, and `--bed`, may be specified
multiple times.

In all cases, 'chr' will be stripped from the chromosome names.

# Interpreting the output

## Interactive example

See: https://brwnj.github.io/covviz/

## Scaled chromosome coverage

Significant regions will be displayed in color atop a gray region which
represents the upper and lower bounds of a given point minus any values
deemed significant.

![significant_regions](data/img/significant_regions.png)

When plotting fewer samples than `--min-samples`, the gray area plot
will not be displayed. Instead, all sample plot traces will be shown.

![min_samples](data/img/min_samples.png)

## Proportions covered

![proportional_coverage](data/img/proportional_coverage.png)

The metadata table will be displayed below the plots.

## Interaction

Clicking on plot traces highlights the line and searches the metadata.
Double-clicking de-selects lines, resets the plot, and de-selects
samples from the table. Clicking on the gene track launches a search
for the gene's respective Gene Card. In cases where genes overlap,
multiple windows/tabs will be opened.

# License

covviz is free and unrestricted for non-commercial use. For commercial use,
please contact [bpedersen@base2genomics.com].
