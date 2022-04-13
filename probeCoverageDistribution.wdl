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
          inputBed = bedFile,
          genomeFile = getGenomeFile.genomeFile,
          outputPrefix = outputFileNamePrefix
      }
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
      coverageHist = calculateProbeCoverageDistribution.coverageHistogram,
      inputBed = bed,
      outputPrefix = outputFileNamePrefix
    }
  }

  output {
    File? coverageHistogram = calculateProbeCoverageDistribution.coverageHistogram
    Array [File]? coverageHistograms = calcProbeCovDistScattered.coverageHistogram
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

task countColumns {
  input {
    File inputBed
    Int jobMemory = 10
    Int timeout = 4
    #String modules = "samtools/1.14"

  }

  parameter_meta {
    inputBed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    jobMemory: "Memory (in GB) allocated for job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    #modules: "Environment module names and version to load (space separated) before command execution."
  }

  meta {
    output_meta : {
      countColumns: "Genome file the defines the expected chromosome order in the bam file."
    }
  }

  command <<<
    cat ~{inputBed} | awk '{print NF}' | sort -nu > number_columns.txt
    #if [[ $numberColumns =   4 ]];
    #then echo "Found a Tomcat!"; fi
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    timeout: "~{timeout}"
    #modules: "~{modules}"
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
    #String modules = "samtools/1.14"

  }

  parameter_meta {
    inputBed: "Target probes, genomic coordinates of the targeted regions in tab-delimited text format."
    jobMemory: "Memory (in GB) allocated for job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    #modules: "Environment module names and version to load (space separated) before command execution."
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
    #modules: "~{modules}"
  }

  output {
    Array[File] splitBeds = glob("*.bed")
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
    Int jobMemory = 10
    Int timeout = 4
    String modules = "rmarkdown/0.1m"

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

  command <<<
    Rscript --vanilla /.mounts/labs/gsiprojects/gsi/gsiusers/blujantoro/wdl/TSprobeCoverage/probeCoverageDistribution/src/plot_coverage_histograms.R \
    -b ~{inputBed} -c  ~{coverageHist} -o ~{outputPrefix}
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
