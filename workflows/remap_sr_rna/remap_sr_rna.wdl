version 1.0

workflow sr_remap_rna {
  input {
    File bam_or_cram
    File bed
    String out_prefix
    Int threads = 1
    File? old_ref
    String star_index
    String genome_dir
  }

  # 1️⃣ Filter reads by BED
  call filter_reads {
    input:
      bam_or_cram = bam_or_cram,
      bed = bed,
      out_prefix = out_prefix,
      threads = threads,
      old_ref = old_ref
  }

  # 2️⃣ Convert filtered BAM to FASTQ
  call extract_reads {
    input:
      bam_or_cram = bam_or_cram,
      bed = bed,
      out_prefix = out_prefix,
      threads = threads, 
      old_ref = old_ref
  }

  # 3️⃣ Remap reads using STAR
  call remap_reads {
    input:
      r1 = extract_reads.r1_fastq,
      r2 = extract_reads.r2_fastq,
      out_prefix = out_prefix,
      threads = threads,
      star_index = star_index,
      genome_dir = genome_dir
  }

  # 4️⃣ Merge filtered (original) and remapped BAMs
  call merge_reads {
    input:
      filtered_bam = filter_reads.filtered_bam,
      remapped_bam = remap_reads.remapped_bam,
      out_prefix = out_prefix,
      threads = threads
  }

  output {
    File final_merged_bam = merge_reads.merged_bam
    File final_merged_bai = merge_reads.merged_bai
  }
}

# -----------------------------------------------------
# Task 1: Filter reads from BAM/CRAM by BED
# -----------------------------------------------------
task filter_reads {
  input {
    File bam_or_cram
    File bed
    String out_prefix
    Int threads
    File? old_ref
    Int memory
    Int disk_size
  }

  command <<<
    
    set -euo pipefail

    echo "Filtering original BAM/CRAM to BED regions..."
    if [[ "~{bam_or_cram}" == *.cram ]]; then
      if [ -z "~{old_ref}" ]; then
        echo "Error: CRAM input requires old reference."
        exit 1
      fi
      samtools view -L "~{bed}" -U "~{out_prefix}_filtered_original.bam" -b --reference "~{old_ref}" -o "~{out_prefix}_filtered_remapped.bam" "~{bam_or_cram}"
      # TODO: Filter transcriptome BAM (?) 
    else
      samtools view -L "~{bed}" -U "~{out_prefix}_filtered_original.bam" -b -o "~{out_prefix}_filtered_remapped.bam" "~{bam_or_cram}" 
      # TODO: Filter transcriptome BAM (?)
    fi

    samtools index "~{out_prefix}_filtered_original.bam"
  >>>

  output {
    File filtered_bam = "~{out_prefix}_filtered_original.bam"
    File filtered_bai = "~{out_prefix}_filtered_original.bam.bai"
  }

  runtime {
    cpu: threads
    memory: "${memory}GB"
    disks: "local-disk ${disk_size} SSD"
    # preemptible: preemptible
    # docker: "biocontainers/samtools:v1.17-4-deb_cv1"
    docker: "staphb/samtools:1.22"
  }
}

# -----------------------------------------------------
# Task 2: Extract reads to FASTQ
# -----------------------------------------------------
task extract_reads {
  input {
    File bam_or_cram
    File bed
    String out_prefix
    Int threads
    File? old_ref
    Int memory
    Int disk_size
  }

  command <<<
    set -euo pipefail

    samtools view -hb -L "~{bed}" -o "~{out_prefix}_extracted_reads.bam" "~{bam_or_cram}"

    echo "Extracting FASTQs from ~{out_prefix}_extracted_reads.bam..."
    samtools fastq \
      -@ "~{threads}" \
      "~{out_prefix}_extracted_reads.bam" \
      -1 "~{out_prefix}.R1.fastq.gz" \
      -2 "~{out_prefix}.R2.fastq.gz" \
      -0 /dev/null -s /dev/null -n
  >>>

  output {
    File r1_fastq = "~{out_prefix}.R1.fastq.gz"
    File r2_fastq = "~{out_prefix}.R2.fastq.gz"
  }

  runtime {
    cpu: threads
    memory: "${memory}GB"
    disks: "local-disk ${disk_size} SSD"
    # preemptible: preemptible
    # docker: "biocontainers/samtools:v1.17-4-deb_cv1"
    docker: "broadgdac/samtools:1.10"
  }
}

# -----------------------------------------------------
# Task 3: Remap reads using STAR
# -----------------------------------------------------
task remap_reads {
  input {
    File r1
    File r2
    String out_prefix
    Int threads
    String star_index
    String genome_dir
    Int memory
    Int disk_size
  }

  command <<<
    set -euo pipefail

    echo "Extracting STAR index"
    mkdir -p "star_index"
    tar -C "star_index" -xzf "~{star_index}" --strip-components=1

    echo "Remapping reads using STAR..."
    STAR \
      --runMode alignReads \
      --runThreadN ~{threads} \
      --readFilesCommand zcat \
      --genomeDir "~{genome_dir}" \
      --readFilesIn "~{r1}" "~{r2}" \
      --outSAMtype "BAM Unsorted" \
      --quantMode "TranscriptomeSAM GeneCounts" \
      # --outFileNamePrefix "~{out_prefix}_remapped_"

    # samtools index "~{out_prefix}_remapped_Aligned.out.bam"

  >>>

  output {
    File remapped_bam = "~{out_prefix}_remapped_Aligned.out.bam"
    # File remapped_bai = "~{out_prefix}_remapped_Aligned.out.bam.bai"
    File remapped_transcriptome_bam = "~{out_prefix}_remapped_Aligned.toTranscriptome.out.bam"
  }

  runtime {
    cpu: threads
    memory: "${memory}GB"
    disks: "local-disk ${disk_size} SSD"
    # docker: "quay.io/biocontainers/star:2.7.10b--h9ee0642_0"
    docker: "us-central1-docker.pkg.dev/depmap-omics/terra-images/star_arriba"
  }
}

# -----------------------------------------------------
# Task 4: Merge original filtered + remapped BAMs
# -----------------------------------------------------
task merge_reads {
  input {
    File filtered_bam
    File remapped_bam
    String out_prefix
    Int threads
    Int memory
    Int disk_size
  }

  command <<<
    set -euo pipefail

    echo "Merging ~{filtered_bam} and ~{remapped_bam}..."
    samtools merge -@ ~{threads} "~{out_prefix}_merged.bam" "~{filtered_bam}" "~{remapped_bam}"
    samtools sort -@ ~{threads} -o "~{out_prefix}_merged.sorted.bam" "~{out_prefix}_merged.bam"
    samtools index "~{out_prefix}_merged.sorted.bam"

    # Add to FixItFelix - TODO: Merge filtered transcriptome BAM to remapped transcriptome BAM
    samtools merge -@ ~{threads} "~{out_prefix}_merged.toTranscriptome.bam" "~{filtered_transcriptome_bam}" "~{remapped_transcriptome_bam}"
  >>>

  output {
    File merged_bam = "~{out_prefix}_merged.sorted.bam"
    File merged_bai = "~{out_prefix}_merged.sorted.bam.bai"
  }

  runtime {
    cpu: threads
    memory: "${memory}GB"
    disks: "local-disk ${disk_size} SSD"
    # docker: "biocontainers/samtools:v1.17-4-deb_cv1"
    # preemptible: preemptible
    # docker: "broadgdac/samtools:1.10"
    docker: "staphb/samtools:1.22"
  }
}
