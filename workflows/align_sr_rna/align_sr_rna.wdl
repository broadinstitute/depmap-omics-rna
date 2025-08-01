version 1.0

workflow align_sr_rna {
    input {
        String sample_id
        String input_file_type
        Array[File]? fastqs
        File? cram_bam
        File? crai_bai
        File? ref_fasta
        File? ref_fasta_index
        Boolean downsample_if_necessary = false
        Int? max_n_reads
        File star_index
    }

    call align_with_star {
        input:
            sample_id = sample_id,
            input_file_type = input_file_type,
            fastqs = fastqs,
            cram_bam = cram_bam,
            crai_bai = crai_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            downsample_if_necessary = downsample_if_necessary,
            max_n_reads = max_n_reads,
            star_index = star_index
    }

    output {
        File analysis_ready_bam = align_with_star.analysis_ready_bam
        File junctions = align_with_star.junctions
        File transcriptome_bam = align_with_star.transcriptome_bam
        File reads_per_gene = align_with_star.reads_per_gene
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
        Boolean downsample_if_necessary = false
        Int? max_n_reads
        File star_index

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/star_arriba"
        String docker_image_hash_or_tag = ":production"
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
                ceil(size(select_first([cram_bam, "/dev/null"]), "GiB") * 8)
            else # BAM
                ceil(size(select_first([cram_bam, "/dev/null"]), "GiB") * 2)
        )
    ) + ceil(size(star_index, "GiB") * 5) + 20 + additional_disk_gb

    Int n_threads = cpu - 1

    command <<<
        set -euo pipefail

        if [[ "~{input_file_type}" == "CRAM" ]]; then
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
        elif [[ "~{input_file_type}" == "BAM" ]]; then
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

        if [[ "~{downsample_if_necessary}" == "true" ]]; then
            if [[ -z "~{max_n_reads}" ]]; then
                echo "ERROR: max_n_reads is required for downsampling" >&2
                exit 1
            fi

            FIRST_FASTQ=$(echo $FASTQS_OPTION | awk '{print $1}')

            if [[ "$FIRST_FASTQ" == *.gz ]]; then
                N_READS=$(zcat "$FIRST_FASTQ" | awk 'END {print NR/4}')
            else
                N_READS=$(awk 'END {print NR/4}' "$FIRST_FASTQ")
            fi

            if (( $(echo "$N_READS > ~{max_n_reads}" | bc -l) )); then
                echo "Downsampling $N_READS reads to ~{max_n_reads}"

                for fq in $FASTQS_OPTION; do
                    if [[ "$fq" == *.gz ]]; then
                        base="${fq%.gz}"
                        zcat "$fq" | seqtk sample -2 -s100 - ~{max_n_reads} \
                            > "${base%.fq}.less.fastq"
                    else
                        seqtk sample -2 -s100 "$fq" ~{max_n_reads} \
                            > "${fq%.fq}.less.fastq"
                    fi

                    mv "${fq%.fq}.less.fastq" "$fq"
                done
            else
                echo "$N_READS <= ~{max_n_reads} reads (no downsampling)"
            fi
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
