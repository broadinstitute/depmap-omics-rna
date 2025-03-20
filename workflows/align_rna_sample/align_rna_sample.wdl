version 1.0

workflow align_rna_sample {
    input {
        String workflow_version = "1.0" # internal semver
        String workflow_source_url # populated automatically with URL of this script

        String sample_id
        String cram_or_bam
        File cram_bam
        File? crai_bai
        File? ref_fasta
        File? ref_fasta_index
        File star_index
    }

    if (cram_or_bam == "BAM") {
        call align_with_star as align_bam_with_star {
            input:
                sample_id = sample_id,
                cram_or_bam = cram_or_bam,
                cram_bam = cram_bam,
                star_index = star_index
        }
    }

    if (cram_or_bam == "CRAM") {
        call align_with_star as align_cram_with_star {
            input:
                sample_id = sample_id,
                cram_or_bam = cram_or_bam,
                cram_bam = cram_bam,
                crai_bai = crai_bai,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                star_index = star_index
        }
    }

    output {
        File transcriptome_bam = select_first([
            align_bam_with_star.transcriptome_bam,
            align_cram_with_star.transcriptome_bam
        ])
        File analysis_ready_bam = select_first([
            align_bam_with_star.analysis_ready_bam,
            align_cram_with_star.analysis_ready_bam
        ])
        File reads_per_gene = select_first([
            align_bam_with_star.reads_per_gene,
            align_cram_with_star.reads_per_gene
        ])
    }
}

task align_with_star {
    input {
        String sample_id
        File cram_bam
        File? crai_bai
        String cram_or_bam
        File? ref_fasta
        File? ref_fasta_index
        File star_index

        String docker_image
        String docker_image_hash_or_tag
        Int cpu = 16
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Float bam_size_gb = if cram_or_bam == "BAM" then (
        size(cram_bam, "GiB")
    ) else 2 * size(cram_bam, "GiB")

    Int disk_space = (
        ceil(bam_size_gb * 10 + size(star_index, "GiB") * 5)
        + 20 + additional_disk_gb
    )

    Int n_threads = cpu - 1

    Int dyn_mem_gb = ceil(10 * bam_size_gb)
    Int max_mem_gb = if (dyn_mem_gb > 128) then 128 else dyn_mem_gb
    Int mem_gb = if (max_mem_gb < 32) then 32 else max_mem_gb

    command <<<
        set -euo pipefail

        echo "Converting CRAM/BAM to FASTQ"
        if [[ "~{cram_or_bam}" == "CRAM" ]];
        then
            samtools sort -n \
                -@ ~{n_threads} \
                --reference "~{ref_fasta}" \
                "~{cram_bam}" \
            | samtools fastq \
                -@ ~{n_threads} \
                --reference "~{ref_fasta}" \
                -1 "~{sample_id}.1.fq.gz" \
                -2 "~{sample_id}.2.fq.gz" -
        else
            samtools sort -n \
                -@ ~{n_threads} \
                "~{cram_bam}" \
            | samtools fastq \
                -@ ~{n_threads} \
                -1 "~{sample_id}.1.fq.gz" \
                -2 "~{sample_id}.2.fq.gz" -
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
            --readFilesIn "~{sample_id}.1.fq.gz" "~{sample_id}.2.fq.gz" \
            --runThreadN ~{n_threads} \
            --twopassMode None \
            --winAnchorMultimapNmax 50
    >>>

    output {
        File analysis_ready_bam = "~{sample_id}.Aligned.out.bam"
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
