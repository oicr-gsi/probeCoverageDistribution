version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem

workflow probeCoverageDistribution {
  input {
    File? fastqR1
    File? fastqR2
    File? bam
    File? bamIndex
    File bed
    String outputFileNamePrefix
    String inputType
    String? partition
  }
  parameter_meta {
    fastqR1: "fastq file for read 1 (optional)."
    fastqR2: "fastq file for read 2 (optional)."
    bam: "Alignment file (optional)."
    bamIndex: "Alignment file index."
    bed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    outputFileNamePrefix: "Optional output prefix to prefix output file names with."
    inputType: "fastq or bam to indicate type of input."
    partition: "Comma separated string indicating pool(s) to be viewed separately in the plots (example: \"exome,pool_1\")."
  }

  if(inputType=="fastq" && defined(fastqR1) && defined(fastqR2)){
    call bwaMem.bwaMem {
      input:
        fastqR1 = select_first([fastqR1]),
        fastqR2 = select_first([fastqR2]),
        outputFileNamePrefix = outputFileNamePrefix,
        readGroups = "'@RG\\tID:ID\\tSM:SAMPLE'",
        doTrim = false
    }
  }

  call getGenomeFile {
    input:
      inputBam = select_first([bwaMem.bwaMemBam,bam])
  }

  #workflow assumes all bed files have 4
  call calculateProbeCoverageDistribution {
    input:
      inputBam = select_first([bwaMem.bwaMemBam,bam]),
      inputBai = select_first([bwaMem.bwaMemIndex,bamIndex]),
      inputBed = bed,
      genomeFile = getGenomeFile.genomeFile,
      outputPrefix = outputFileNamePrefix
  }


  if(!defined(partition)){
    call Rplot {
      input:
        coverageHist = calculateProbeCoverageDistribution.coverageHistogram,
        inputBed = bed,
        outputPrefix = outputFileNamePrefix
    }
  }

  if(defined(partition)){
    call Rplot as RplotPartioned {
      input:
        coverageHist = calculateProbeCoverageDistribution.coverageHistogram,
        inputBed = bed,
        outputPrefix = outputFileNamePrefix,
        partition = partition
    }
  }

  call compressResults{
    input:
      inFiles= select_first([Rplot.Rplots,RplotPartioned.Rplots]),
      outputPrefix = outputFileNamePrefix
  }

  output {
    File cvgFile = calculateProbeCoverageDistribution.coverageHistogram
    File plots = compressResults.plots
  }

  meta {
    author: "Beatriz Lujan Toro"
    email: "beatriz.lujantoro@oicr.on.ca"
    description: "Workflow to calculate probe coverage"
    dependencies: [{
      name: "bedtools/2.27",
      url: "https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html"
    },
    {
      name: "samtools/1.14",
      url: "http://www.htslib.org/"
    },
    {
      name: "optparse",
      url: "https://cran.r-project.org/web/packages/optparse/index.html"
    },
    {
      name: "ggpubr",
      url: "https://cran.r-project.org/web/packages/ggpubr/index.html"
    },
    {
      name: "ggplot2",
      url: "https://cran.r-project.org/web/packages/ggplot2/index.html"
    },
    {
      name: "tidytext",
      url: "https://cran.r-project.org/web/packages/tidytext/"
    },
    {
      name: "r/3.6.1",
      url: "https://cran.r-project.org/"
    }]
    output_meta: {
      cvgFile: "Coverage histogram, tab-delimited text file reporting the coverage at each feature in the bed file.",
      plots: "A compress file of all the Rplots created by the workflow, which show interval panel coverage."
    }
  }
}

task getGenomeFile {
  input {
    File inputBam
    Int jobMemory = 10
    Int timeout = 4
    String modules = "samtools/1.14"

  }

  parameter_meta {
    inputBam: "Input file (bam)."
    jobMemory: "Memory (in GB) allocated for job."
    modules: "Environment module names and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta : {
      genomeFile: "Genome file the defines the expected chromosome order in the bam file."
    }
  }

  command <<<
    samtools view -H ~{inputBam} | \
    grep @SQ | sed 's/@SQ\tSN:\|LN://g' \
    > genome.txt
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File genomeFile = "genome.txt"
  }
}

task calculateProbeCoverageDistribution {
  input {
    File inputBam
    File inputBai
    File inputBed
    File genomeFile
    Int jobMemory = 10
    Int timeout = 4
    String outputPrefix
    String modules = "bedtools/2.27"
  }

  parameter_meta {
    inputBam: "Input file (bam)."
    inputBai: "index of the input .bam file"
    inputBed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    genomeFile: "Genome file the defines the expected chromosome order in the bam file."
    jobMemory: "Memory (in GB) allocated for job."
    outputPrefix: "Output prefix to prefix output file names with."
    modules: "Environment module names and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta : {
      coverageHistogram: "Coverage histogram, tab-delimited text file reporting the coverage at each feature in the bed file."
    }
  }


  command <<<
    bedtools coverage -hist \
    -a ~{inputBed} \
    -b ~{inputBam} \
    -sorted -g ~{genomeFile} \
    > ~{outputPrefix}.cvghist.txt || echo "Bedtools failed to produce output" \
    | rm ~{outputPrefix}.cvghist.txt
    #use or "||" when command fails otherwise workflow succeeds on empty file
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File coverageHistogram = "~{outputPrefix}.cvghist.txt"
  }
}

task Rplot {
  input {
    File inputBed
    File coverageHist
    String outputPrefix
    String? partition
    Int jobMemory = 20
    Int timeout = 4
    String modules = "probe-coverage-distribution/2.0"
  }

  parameter_meta {
    inputBed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    coverageHist: "Coverage histogram, tab-delimited text file reporting the coverage at each feature in the bed file."
    outputPrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) allocated for job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module names and version to load (space separated) before command execution."
    partition: "Comma separated string indicating pool(s) to be viewed separately in the plots (example: \"exome,pool_1\")."
  }

  meta {
    output_meta : {
      plots: "Probe coverage distribution plots."
    }
  }

  command <<<
    Rscript --vanilla $PROBE_COVERAGE_DISTRIBUTION_ROOT/plot_coverage_histograms.R \
    -b ~{inputBed} -c ~{coverageHist} -o ~{outputPrefix} ~{"-l " + partition}
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    Array[File] Rplots = glob("*.png")
  }
}

task compressResults {

  input {
    Array[File] inFiles
    String outputPrefix
    Int jobMemory = 12
    Int timeout = 4
  }

  parameter_meta {
    inFiles: "Array of input files"
    jobMemory: "Memory for the task, in gigabytes"
    timeout: "Timeout for the task, in hours"
  }

  meta {
    description: "Gather plots into a .tar.gz"
    output_meta: {
      plots: "Compressed folder with plots"
    }
  }

  # create a directory to compress
  command <<<
    set -euo pipefail
    mkdir ~{outputPrefix}
    cp -t ~{outputPrefix} ~{sep=' ' inFiles}
    tar -zcvf "~{outputPrefix}_plots.tar.gz" "~{outputPrefix}"
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    timeout: "~{timeout}"
  }

  output {
    File plots = "~{outputPrefix}_plots.tar.gz"
  }

}
