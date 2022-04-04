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
  }
  parameter_meta {
    fastqR1: "fastq file for read 1"
    fastqR2: "fastq file for read 2"
    bed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    outputFileNamePrefix: "Optional output prefix to prefix output file names with."
  }

  if(inputType=="fastq" && defined(fastqR1) && defined(fastqR2)){
    call bwaMem.bwaMem {
      input:
        fastqR1 = select_first([fastqR1]),
        fastqR2 = select_first([fastqR2]),
        outputFileNamePrefix = outputFileNamePrefix,
        readGroups = "'@RG\\tID:ID\\tSM:SAMPLE'",
        doTrim = false #TEST CHECK LATER
    }
  }

  call getGenomeFile {
    input:
      inputBam = select_first([bwaMem.bwaMemBam,bam])
  }

  call calculateProbeCoverageDistribution {
    input:
      inputBam = select_first([bwaMem.bwaMemBam,bam]),
      inputBai = select_first([bwaMem.bwaMemIndex,bamIndex]),
      inputBed = bed,
      genomeFile = getGenomeFile.genomeFile,
      outputPrefix = outputFileNamePrefix
  }

  output {
    File coverageHistogram = calculateProbeCoverageDistribution.coverageHistogram
  }

  meta {
    author: "Beatriz Lujan Toro"
    email: "beatriz.lujantoro@oicr.on.ca"
    description: "Workflow to calculate probe coverage"
    dependencies: [{
      name: "bedtools/2.27",
      url: "https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html"
    }]
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
    > "~{outputPrefix}.cvghist.txt" || echo "Bedtools failed to produce output" \
    | rm "~{outputPrefix}.cvghist.txt"
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
