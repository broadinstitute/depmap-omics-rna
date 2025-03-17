version 1.0

workflow align_rna_sample {
    input {
        String workflow_version = "1.0" # internal semver
        String workflow_source_url # populated automatically with URL of this script

        String sample_id
        String cram_or_bam
        File cram_bam
        File star_index
    }

    call align_with_star {
        input:
            sample_id = sample_id,
            cram_or_bam = cram_or_bam,
            cram_bam = cram_bam,
            star_index = star_index
    }

    output {
        File transcriptome_bam = align_with_star.transcriptome_bam
        File analysis_ready_bam = align_with_star.analysis_ready_bam
        File analysis_ready_bai = align_with_star.analysis_ready_bai
        File reads_per_gene = align_with_star.reads_per_gene
    }
}

task align_with_star {
    input {
        String sample_id
        File cram_bam
        String cram_or_bam
        File star_index

        String docker_image
        String docker_image_hash_or_tag
        Int cpu = 8
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
        Int additional_mem_gb = 0
    }

    Float bam_size_gb = if cram_or_bam == "BAM" then (
        size(cram_bam, "GiB")
    ) else 2 * size(cram_bam, "GiB")

    Int disk_space = (
        ceil(bam_size_gb * 10 + size(star_index, "GiB") * 5)
        + 20 + additional_disk_gb
    )

    Int n_threads = cpu - 1

    Int mem_gb = ceil(20 * bam_size_gb) + additional_mem_gb
    Int bam_sort_ram_bytes = ceil(mem_gb * 1000 * 1000 * 1000 * 0.85)

    command <<<
        set -euo pipefail

        echo "Extracting STAR index"
        mkdir -p "star_index"
        tar -C "star_index" -xzf "~{star_index}" --strip-components=1

        echo "Converting CRAM/BAM to FASTQ"
        samtools sort -n \
            -@ ~{n_threads} \
            "~{cram_bam}" \
        | samtools fastq \
            -@ ~{n_threads} \
            -1 "~{sample_id}.1.fq.gz" \
            -2 "~{sample_id}.2.fq.gz" -

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
            --limitBAMsortRAM ~{bam_sort_ram_bytes} \
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

        echo "Sorting and indexing aligned BAM"
        samtools sort \
            -@ ~{n_threads} \
            "~{sample_id}.Aligned.out.bam" \
            > "~{sample_id}.analysis_ready.bam"

        samtools index \
            -b \
            -@ ~{n_threads} \
            "~{sample_id}.analysis_ready.bam" \
            "~{sample_id}.analysis_ready.bai"
    >>>

    output {
        File analysis_ready_bai = "~{sample_id}.analysis_ready.bai"
        File analysis_ready_bam = "~{sample_id}.analysis_ready.bam"
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
