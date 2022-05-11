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

  call countColumns {
    input:
      inputBed = bed
  }

  if (countColumns.numberColumns == 5) {
    #workflow assumes all bed files have 4 or 5 columns
    call splitBed {
      input:
        inputBed = bed
    }

    scatter (bedFile in splitBed.splitBeds) {

      call calculateProbeCoverageDistribution as calcProbeCovDistScattered {

        input:
          inputBam = select_first([bwaMem.bwaMemBam,bam]),
          inputBai = select_first([bwaMem.bwaMemIndex,bamIndex]),
          inputBed = bedFile.right,
          genomeFile = getGenomeFile.genomeFile,
          outputPrefix = outputFileNamePrefix
      }

      call Rplot as RplotScattered {
        input:
          coverageHist = calcProbeCovDistScattered.coverageHistogram,
          inputBed = bedFile.right,
          outputPrefix = outputFileNamePrefix,
          multipleBed = bedFile.left
      }

    }

    call zipResults as zipScatteredResults {
     input:
       inFiles=flatten(RplotScattered.Rplots),
       outputPrefix = outputFileNamePrefix
    }
  }

  if (countColumns.numberColumns == 4) {
    #workflow assumes all bed files have 4 or 5 columns
    call calculateProbeCoverageDistribution {
      input:
        inputBam = select_first([bwaMem.bwaMemBam,bam]),
        inputBai = select_first([bwaMem.bwaMemIndex,bamIndex]),
        inputBed = bed,
        genomeFile = getGenomeFile.genomeFile,
        outputPrefix = outputFileNamePrefix
    }

    call Rplot {
      input:
        coverageHist = calculateProbeCoverageDistribution.coverageHistogram,
        inputBed = bed,
        outputPrefix = outputFileNamePrefix,
        #multipleBed = false
    }

    call zipResults{
      input: inFiles=Rplot.Rplots,
      outputPrefix = outputFileNamePrefix
    }
  }

  output {
    Array[File]? cvgFiles = calcProbeCovDistScattered.coverageHistogram
    File? cvgFile = calculateProbeCoverageDistribution.coverageHistogram
    File plots = select_first([zipScatteredResults.zipArchive, zipResults.zipArchive])
  }

  meta {
    author: "Beatriz Lujan Toro"
    email: "beatriz.lujantoro@oicr.on.ca"
    description: "Workflow to calculate probe coverage"
    dependencies: [{
      name: "bedtools/2.27",
      url: "https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html"
    }]
    output_meta: {
      cvgFile: "Coverage histogram, tab-delimited text file reporting the coverage at each feature in the bed file.",
      plots: "A zip file of all the Rplots created by the workflow, which show interval panel coverage."
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

task countColumns {
  input {
    File inputBed
    Int jobMemory = 10
    Int timeout = 4
  }

  parameter_meta {
    inputBed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    jobMemory: "Memory (in GB) allocated for job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta : {
      countColumns: "Genome file the defines the expected chromosome order in the bam file."
    }
  }

  command <<<
    cat ~{inputBed} | awk '{print NF}' | sort -nu > number_columns.txt
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    timeout: "~{timeout}"
  }

  output {
    Int numberColumns = read_int("number_columns.txt")
  }
}

task splitBed {
  input {
    File inputBed
    Int jobMemory = 10
    Int timeout = 4
  }

  parameter_meta {
    inputBed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    jobMemory: "Memory (in GB) allocated for job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta : {
      splitBeds: "Bed file split when there are two name columns (feature description)."
    }
  }

  command <<<
    cat ~{inputBed} | cut -f 1-4 > probes_1.bed
    cat ~{inputBed} | cut -f 1-3,5 > probes_2.bed
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    timeout: "~{timeout}"
  }

  output {
    Array[Pair[String, File]] splitBeds = [("1",  "probes_1.bed"),("2",  "probes_2.bed")]
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

task Rplot {
  input {
    File inputBed
    File coverageHist
    String outputPrefix
    String? multipleBed
    Int jobMemory = 20
    Int timeout = 4
    String modules = "probe-coverage-distribution/1.0"
  }

  parameter_meta {
    inputBed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    coverageHist: "Coverage histogram, tab-delimited text file reporting the coverage at each feature in the bed file."
    outputPrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) allocated for job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module names and version to load (space separated) before command execution."
  }

  meta {
    output_meta : {
      plots: "Probe coverage distribution plots."
    }
  }

  #String s = "~{if b then '${1 + i}' else 0}"

  command <<<
    Rscript --vanilla /.mounts/labs/gsiprojects/gsi/gsiusers/blujantoro/wdl/TSprobeCoverage/probeCoverageDistribution/src/plot_coverage_histograms.R \
    -b ~{inputBed} -c ~{coverageHist} -o ~{outputPrefix}~{"_" + multipleBed}
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

task zipResults {

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
    description: "Gather plots into a .zip archive"
    output_meta: {
      zipArchive: "ZIP archive file"
    }
  }

  # create a directory for the zip archive; allows unzip without exploding multiple files into the working directory

  command <<<
    set -euo pipefail
    mkdir ~{outputPrefix}
    cp -t ~{outputPrefix} ~{sep=' ' inFiles}
    zip -qr ~{outputPrefix}.zip ~{outputPrefix}
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    timeout: "~{timeout}"
  }

  output {
    File zipArchive = "~{outputPrefix}.zip"
  }

}
