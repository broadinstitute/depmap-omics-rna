version 1.0

workflow align_rna_sample {
    input {
        String sample_id
        String input_file_type
        Array[File]? fastqs
        File? cram_bam
        File? crai_bai
        File? ref_fasta
        File? ref_fasta_index
        File star_index
    }

    if (input_file_type == "FASTQ") {
        call align_with_star as align_fastq_with_star {
            input:
                sample_id = sample_id,
                input_file_type = input_file_type,
                fastqs = fastqs,
                star_index = star_index
        }
    }

    if (input_file_type == "BAM") {
        call align_with_star as align_bam_with_star {
            input:
                sample_id = sample_id,
                input_file_type = input_file_type,
                cram_bam = cram_bam,
                star_index = star_index
        }
    }

    if (input_file_type == "CRAM") {
        call align_with_star as align_cram_with_star {
            input:
                sample_id = sample_id,
                input_file_type = input_file_type,
                cram_bam = cram_bam,
                crai_bai = crai_bai,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                star_index = star_index
        }
    }

    output {
        File analysis_ready_bam = select_first([
            align_fastq_with_star.analysis_ready_bam,
            align_bam_with_star.analysis_ready_bam,
            align_cram_with_star.analysis_ready_bam
        ])
        File junctions = select_first([
            align_fastq_with_star.junctions,
            align_bam_with_star.junctions,
            align_cram_with_star.junctions
        ])
        File transcriptome_bam = select_first([
            align_fastq_with_star.transcriptome_bam,
            align_bam_with_star.transcriptome_bam,
            align_cram_with_star.transcriptome_bam
        ])
        File reads_per_gene = select_first([
            align_fastq_with_star.reads_per_gene,
            align_bam_with_star.reads_per_gene,
            align_cram_with_star.reads_per_gene
        ])
    }
}

task align_with_star {
    input {
        String sample_id
        String input_file_type
        Array[File]? fastqs
        File? cram_bam
        File? crai_bai
        File? ref_fasta
        File? ref_fasta_index
        File star_index

        String docker_image
        String docker_image_hash_or_tag
        Int cpu = 16
        Int mem_gb = 48
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        if input_file_type == "FASTQ" then
            ceil(size(select_first([fastqs, []]), "GiB"))
        else (
            if input_file_type == "CRAM" then
                ceil(size(select_first([cram_bam, "/dev/null"]), "GiB") * 3)
            else # BAM
                ceil(size(select_first([cram_bam, "/dev/null"]), "GiB"))
        )
    ) + ceil(size(star_index, "GiB") * 5) + 20 + additional_disk_gb

    Int n_threads = cpu - 1

    command <<<
        set -euo pipefail

        if [[ "~{input_file_type}" == "CRAM" ]];
        then
            echo "Converting CRAM to FASTQ"

            samtools sort -n \
                -@ ~{n_threads} \
                --reference "~{ref_fasta}" \
                "~{cram_bam}" \
            | samtools fastq \
                -@ ~{n_threads} \
                --reference "~{ref_fasta}" \
                -1 "~{sample_id}.1.fq.gz" \
                -2 "~{sample_id}.2.fq.gz" -

            FASTQS_OPTION="~{sample_id}.1.fq.gz ~{sample_id}.2.fq.gz"
        elif [[ "~{input_file_type}" == "BAM" ]];
        then
            echo "Converting BAM to FASTQ"

            samtools sort -n \
                -@ ~{n_threads} \
                "~{cram_bam}" \
            | samtools fastq \
                -@ ~{n_threads} \
                -1 "~{sample_id}.1.fq.gz" \
                -2 "~{sample_id}.2.fq.gz" -

            FASTQS_OPTION="~{sample_id}.1.fq.gz ~{sample_id}.2.fq.gz"
        else
            FASTQS_OPTION="~{sep=' ' fastqs}"
        fi

        echo "Extracting STAR index"
        mkdir -p "star_index"
        tar -C "star_index" -xzf "~{star_index}" --strip-components=1

        echo "Running STAR"
        STAR \
            --alignInsertionFlush Right \
            --alignIntronMax 100000 \
            --alignMatesGapMax 100000 \
            --alignSJDBoverhangMin 10 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --alignSplicedMateMapLmin 30 \
            --alignSplicedMateMapLminOverLmate 0 \
            --chimJunctionOverhangMin 8 \
            --chimMultimapNmax 20 \
            --chimMultimapScoreRange 3 \
            --chimNonchimScoreDropMin 10 \
            --chimOutJunctionFormat 1 \
            --chimOutType Junctions WithinBAM \
            --chimScoreJunctionNonGTAG -4 \
            --chimSegmentMin 12 \
            --genomeDir star_index \
            --genomeLoad NoSharedMemory \
            --limitSjdbInsertNsj 1200000 \
            --outFileNamePrefix "~{sample_id}." \
            --outFilterIntronMotifs None \
            --outFilterMatchNmin 0 \
            --outFilterMatchNminOverLread 0.33 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.1 \
            --outFilterMultimapNmax 20 \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterType BySJout \
            --outReadsUnmapped None \
            --outSAMattrRGline ID:GRPundef \
            --outSAMattributes NH HI AS nM NM ch \
            --outSAMstrandField intronMotif \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --peOverlapMMp 0.1 \
            --peOverlapNbasesMin 12 \
            --quantMode TranscriptomeSAM GeneCounts \
            --readFilesCommand gunzip -c \
            --readFilesIn $FASTQS_OPTION \
            --runThreadN ~{n_threads} \
            --twopassMode None \
            --winAnchorMultimapNmax 50
    >>>

    output {
        File analysis_ready_bam = "~{sample_id}.Aligned.out.bam"
        File junctions = "~{sample_id}.SJ.out.tab"
        File reads_per_gene = "~{sample_id}.ReadsPerGene.out.tab"
        File transcriptome_bam = "~{sample_id}.Aligned.toTranscriptome.out.bam"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    meta {
        allowNestedInputs: true
    }
}
