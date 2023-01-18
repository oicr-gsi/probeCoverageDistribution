# probeCoverageDistribution

Workflow to visualize probe coverage

## Overview

## Dependencies

* [bedtools 2.27](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
* [samtools 1.14](http://www.htslib.org/)
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
* [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [tidytext](https://cran.r-project.org/web/packages/tidytext/)
* [r 3.6.1](https://cran.r-project.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run probeCoverageDistribution.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`bed`|File|Target probes, genomic coordinates of the targeted regions in tab-delimited text format.
`outputFileNamePrefix`|String|Optional output prefix to prefix output file names with.
`inputType`|String|fastq or bam to indicate type of input.
`bwaMem.runBwaMem_bwaRef`|String|The reference genome to align the sample with by BWA
`bwaMem.runBwaMem_modules`|String|Required environment modules


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR1`|File?|None|fastq file for read 1 (optional).
`fastqR2`|File?|None|fastq file for read 2 (optional).
`bam`|File?|None|Alignment file (optional).
`bamIndex`|File?|None|Alignment file index.
`partition`|String?|None|Comma separated string indicating pool(s) to be viewed separately in the plots (example: "exome,pool_1").


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`bwaMem.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.indexBam_timeout`|Int|48|Hours before task timeout
`bwaMem.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwaMem.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.bamMerge_timeout`|Int|72|Hours before task timeout
`bwaMem.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwaMem.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwaMem.runBwaMem_timeout`|Int|96|Hours before task timeout
`bwaMem.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`bwaMem.runBwaMem_threads`|Int|8|Requested CPU threads
`bwaMem.runBwaMem_addParam`|String?|None|Additional BWA parameters
`bwaMem.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`getGenomeFile.jobMemory`|Int|10|Memory (in GB) allocated for job.
`getGenomeFile.timeout`|Int|4|Maximum amount of time (in hours) the task can run for.
`getGenomeFile.modules`|String|"samtools/1.14"|Environment module names and version to load (space separated) before command execution.
`calculateProbeCoverageDistribution.jobMemory`|Int|10|Memory (in GB) allocated for job.
`calculateProbeCoverageDistribution.timeout`|Int|4|Maximum amount of time (in hours) the task can run for.
`calculateProbeCoverageDistribution.modules`|String|"bedtools/2.27"|Environment module names and version to load (space separated) before command execution.
`Rplot.partition`|String?|None|Comma separated string indicating pool(s) to be viewed separately in the plots (example: "exome,pool_1").
`Rplot.jobMemory`|Int|20|Memory (in GB) allocated for job.
`Rplot.timeout`|Int|4|Maximum amount of time (in hours) the task can run for.
`Rplot.modules`|String|"probe-coverage-distribution/2.0"|Environment module names and version to load (space separated) before command execution.
`RplotPartioned.jobMemory`|Int|20|Memory (in GB) allocated for job.
`RplotPartioned.timeout`|Int|4|Maximum amount of time (in hours) the task can run for.
`RplotPartioned.modules`|String|"probe-coverage-distribution/2.0"|Environment module names and version to load (space separated) before command execution.
`compressResults.jobMemory`|Int|12|Memory for the task, in gigabytes
`compressResults.timeout`|Int|4|Timeout for the task, in hours


### Outputs

Output | Type | Description
---|---|---
`cvgFile`|File|Coverage histogram, tab-delimited text file reporting the coverage at each feature in the bed file.
`plots`|File|A compress file of all the Rplots created by the workflow, which show interval panel coverage.


## Commands
 This section lists command(s) run by the probeCoverageDistribution workflow
 
 * Running probeCoverageDistribution workflow
 
 Capture panels are defined by a set of intervals defined in a bed file. This workflow uses bedtools to obtain bam coverage of the set of intervals and produces summary plots to visualize probe coverage and assess performance. 
 
 Create genome file to use bedtools sorted to estimate coverage
 ```
     samtools view -H ~{inputBam} | \
     grep @SQ | sed 's/@SQ\tSN:\|LN://g' \
     > genome.txt
 ```
 Calculate coverage using bedtools
 ```
     bedtools coverage -hist \
     -a ~{inputBed} \
     -b ~{inputBam} \
     -sorted -g ~{genomeFile} \
     > ~{outputPrefix}.cvghist.txt || echo "Bedtools failed to produce output" \
     | rm ~{outputPrefix}.cvghist.txt
     #use or "||" when command fails otherwise workflow succeeds on empty file
 ```
 Run custom Rscript to plot coverage
 ```
     Rscript --vanilla $PROBE_COVERAGE_DISTRIBUTION_ROOT/plot_coverage_histograms.R \
     -b ~{inputBed} -c ~{coverageHist} -o ~{outputPrefix} ~{"-l " + partition}
 ```
 Compress plots
 ```
     set -euo pipefail
     mkdir ~{outputPrefix}
     cp -t ~{outputPrefix} ~{sep=' ' inFiles}
     tar -zcvf "~{outputPrefix}_plots.tar.gz" "~{outputPrefix}"
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
