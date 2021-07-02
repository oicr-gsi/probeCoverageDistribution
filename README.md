# probeCoverageDistribution

Workflow to calculate probe coverage

## Overview

## Dependencies

* [bedtools 2.27](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)


## Usage

### Cromwell
```
java -jar cromwell.jar run probeCoverageDistribution.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputBed`|File|Target probes, genomic coordinates of the targeted regions in tab-delimited text format.
`inputBam`|File|Input file bam.


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|basename(inputBam,'.bam')|Output prefix to prefix output file names with.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`calculateProbeCoverageDistribution.jobMemory`|Int|10|Memory (in GB) allocated for job.
`calculateProbeCoverageDistribution.timeout`|Int|4|Maximum amount of time (in hours) the task can run for.
`calculateProbeCoverageDistribution.modules`|String|"bedtools/2.27"|Environment module names and version to load (space separated) before command execution.


### Outputs

Output | Type | Description
---|---|---
`coverageHistogram`|File|Coverage histogram, tab-delimited text file reporting the coverage at each feature in the bed file.


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
