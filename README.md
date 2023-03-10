<!-- # ![nf-core/manticore](docs/images/nf-core-manticore_logo_light.png#gh-light-mode-only) ![nf-core/manticore](docs/images/nf-core-manticore_logo_dark.png#gh-dark-mode-only) -->

<!-- [![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/manticore/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<!-- [![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/manticore) -->

<!-- [![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23manticore-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/manticore)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core) -->

## Introduction

**nf-core/manticore** is a bioinformatics best-practice analysis
pipeline for Population genomics analyses of non-model organisms.

The pipeline is built using [Nextflow](https://www.nextflow.io), a
workflow tool to run tasks across multiple compute infrastructures in
a very portable manner. It uses Docker/Singularity containers making
installation trivial and results highly reproducible. The [Nextflow
DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of
this pipeline uses one container per process which makes it much
easier to maintain and update software dependencies. Where possible,
these processes have been submitted to and installed from
[nf-core/modules](https://github.com/nf-core/modules) in order to make
them available to all nf-core pipelines, and to everyone within the
Nextflow community!

## Pipeline summary

1. Calculate coverage distribution from bam files with [`mosdepth`](https://github.com/brentp/mosdepth)
2. Generate sequence masks from coverage information
3. WIP: Calculate genetic variation (pi, theta, S) with [`vcftools`](https://vcftools.sourceforge.net/) for different sequence masks and window sizes.
4. WIP: Calculate within and between population differentiation statistics (FST) with [`vcftools`](https://vcftools.sourceforge.net/) for user-defined sample set configuration
5. TODO: Make a preliminary principal component analysis (PCA) with [`plink`](https://www.cog-genomics.org/plink/2.0/) to assess population structure
6. TODO: Estimate admixture components with [`ADMIXTURE`](https://dalexander.github.io/admixture/index.html)
7. WIP: Run windowed selection scans (Tajima's D)

## Input

The workflow requires two inputs:

1. vcf file with variants
2. samplesheet - a comma-separated file with columns `sample_id`, `bam` file path, `bai` file path, e.g.,

```bash
sample,bam,bai
YRI-1,data/YRI-1.chr22.bam,data/YRI-1.chr22.bam.bai
YRI-2,data/YRI-2.chr22.bam,data/YRI-2.chr22.bam.bai
CEU-1,data/CEU-1.chr22.bam,data/CEU-1.chr22.bam.bai
CEU-2,data/CEU-2.chr22.bam,data/CEU-2.chr22.bam.bai
```

### Sample sets

Samples can be grouped into sample sets (aka populations) and passed
to the `--sample_sets` option. The sample set file is a tab-separated
file consisting of two columns `sampleset_id` and a comma-separated
column of sample names, e.g.,

    YRI    YRI-1,YRI-2
    CEU    CEU-1,CEU-2

Sample sets listed in the input file will be pairwise compared where
applicable (e.g., fst).

A default sample set `ALL`, consisting of all samples, is
automatically constructed.

### Coverage filters

Coverage filters are subsequently applied
individually to all sample sets. An automatic coverage filter, denoted
`auto`, is automatically applied where positions having coverage
within the range coverage mean +/- 0.5 standard deviations are saved
to bed files denoting regions that pass filtering. In addition, the
options `--coverage_min` and `--coverage_max` allow setting manual
filters, denoted `manual`, on a sampleset specific basis, e.g.,

    --coverage_min ALL:3,YRI:2 --coverage_max ALL:20

### Window analyses

The option `--window_sizes` takes a comma-separated string of integers
corresponding to window sizes what will be applied to window-based
analyses:

    --window_sizes 100000,500000

### Regions of interest (ROI)

The option `--roi_fof` is a comma-separated file consisting of a
column denoting the type of analysis (`site` or `window`) and a path
to a bed file defining regions of interest:

    site,tests/resources/exons.bed
    window,tests/resources/utr.bed

Results are calculated for the ROI intersected with the coverage
filtered regions.

A default ROI, `genome`, corresponding to the entire genome, is added
automatically for whole-genome analyses.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

```bash
nextflow run nf-core/manticore -profile test,YOURPROFILE --outdir <OUTDIR>
```

Note that some form of configuration will be needed so that
Nextflow knows how to fetch the required software. This is usually
done in the form of a config profile (`YOURPROFILE` in the example
command above). You can chain multiple config profiles in a
comma-separated string.

> - The pipeline comes with config profiles called `docker`,
>   `singularity`, `podman`, `shifter`, `charliecloud` and `conda`
>   which instruct the pipeline to use the named tool for software
>   management. For example, `-profile test,docker`. - Please check
>   [nf-core/configs](https://github.com/nf-core/configs#documentation)
>   to see if a custom config file to run nf-core pipelines already
>   exists for your Institute. If so, you can simply use `-profile
<institute>` in your command. This will enable either `docker`
>   or `singularity` and set the appropriate execution settings for
>   your local compute environment. - If you are using
>   `singularity`, please use the [`nf-core
download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use)
>   command to download images first, before running the pipeline.
>   Setting the [`NXF_SINGULARITY_CACHEDIR` or
>   `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub)
>   Nextflow options enables you to store and re-use the images
>   from a central location for future pipeline runs. - If you are
>   using `conda`, it is highly recommended to use the
>   [`NXF_CONDA_CACHEDIR` or
>   `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html)
>   settings to store the environments in a central location for
>   future pipeline runs.

4.  Start running your own analysis!

```bash
nextflow run nf-core/manticore --input samplesheet.csv --outdir <OUTDIR> --fasta reference.fasta -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
```

## Documentation

<!-- The nf-core/manticore pipeline comes with documentation about the pipeline [usage](https://nf-co.re/manticore/usage), [parameters](https://nf-co.re/manticore/parameters) and [output](https://nf-co.re/manticore/output). -->

## Credits

nf-core/manticore was originally written by Per Unneberg.

<!-- We thank the following people for their extensive assistance in the -->
<!-- development of this pipeline: -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the
[contributing guidelines](.github/CONTRIBUTING.md).

<!-- For further information or help, don't hesitate to get in touch on the -->
<!-- [Slack `#manticore` -->
<!-- channel](https://nfcore.slack.com/channels/manticore) (you can join -->
<!-- with [this invite](https://nf-co.re/join/slack)). -->

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/manticore for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
