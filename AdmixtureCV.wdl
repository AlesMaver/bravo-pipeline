version 1.0
## Copyright CMG@KIGM, Peter Juvan

# WORKFLOW DEFINITION 
workflow AdmixtureCV {
  input {
    File input_bed
    File input_bim
    File input_fam
    Array[Int] num_populations = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    Int threads = 256
    Int cv_fold = 5
    String output_basename = basename(input_bed, ".bed")
  }

  scatter (K in num_populations) {

    call crossValAdmixture {
      input:
        input_bed = input_bed,
        input_bim = input_bim,
        input_fam = input_fam,
        threads = threads,
        K = K,
        cv_fold = cv_fold,
        output_basename = output_basename
    }

  } # Close per K scatter

  call collectCVerror {
    input:
      inputs_stdout = crossValAdmixture.output_stdout,
      output_basename = output_basename
  }

  output {
    Array[File] output_P = crossValAdmixture.output_P
    Array[File] output_Q = crossValAdmixture.output_Q
    File output_CVerror = collectCVerror.output_CVerror
  }
} # Close workflow


task crossValAdmixture {
  input {
    File input_bed
    File input_bim
    File input_fam
    Int threads
    Int K
    Int cv_fold
    String output_basename
  }

  command {
    set -e
    admixture --cv=~{cv_fold} -j~{threads} ~{input_bed} ~{K} | tee "~{output_basename}.~{K}.stdout"
  }

  runtime {
    docker: "quay.io/biocontainers/admixture:1.3.0--0"
    #requested_memory_mb_per_core: 256
    cpu: threads
  }
  output {
    File output_P = "~{output_basename}.~{K}.P"
    File output_Q = "~{output_basename}.~{K}.Q"
    File output_stdout = "~{output_basename}.~{K}.stdout"
  }
}

task collectCVerror {
  input {
    Array[File] inputs_stdout
    String output_basename
  }

  command <<<
    set -e
    grep -h CV ~{sep=' ' inputs_stdout} | mawk '{gsub("\(K\=", ""); gsub("\)\:", ""); print($3,$4)}' | tee "~{output_basename}.CVerror.txt"
  >>>

  runtime {
    docker: "ubuntu"
    requested_memory_mb_per_core: 1000
    cpu: 1
  }
  output {
    File output_CVerror = "~{output_basename}.CVerror.txt"
  }
}