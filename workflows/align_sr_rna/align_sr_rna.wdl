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
        File star_index
        String align_insertion_flush = "Right"
        Int align_intron_max = 100000
        Int align_mates_gap_max = 100000
        Int align_sjdb_overhang_min = 10
        String align_sj_stitch_mismatch_nmax = "5 -1 5 5"
        Int align_spliced_mate_map_lmin = 30
        Int align_spliced_mate_map_lmin_over_lmate = 0
        Int chim_junction_overhang_min = 8
        Int chim_multimap_nmax = 20
        Int chim_multimap_score_range = 3
        Int chim_nonchim_score_drop_min = 10
        Int chim_out_junction_format = 1
        String chim_out_type = "Junctions WithinBAM"
        Int chim_score_junction_non_gtag = -4
        Int chim_segment_min = 12
        String genome_dir = "star_index"
        String genome_load = "NoSharedMemory"
        Int limit_sjdb_insert_nsj = 1200000
        String out_filter_intron_motifs = "None"
        Int out_filter_match_nmin = 0
        Float out_filter_match_nmin_over_lread = 0.33
        Int out_filter_mismatch_nmax = 999
        Float out_filter_mismatch_nover_lmax = 0.1
        Int out_filter_multimap_nmax = 20
        Float out_filter_score_min_over_lread = 0.33
        String out_filter_type = "BySJout"
        String out_reads_unmapped = "None"
        String out_sam_attr_rg_line = "ID:GRPundef"
        String out_sam_attributes = "NH HI AS nM NM ch"
        String out_sam_strand_field = "intronMotif"
        String out_sam_type = "BAM Unsorted"
        String out_sam_unmapped = "Within"
        Float pe_overlap_mmp = 0.1
        Int pe_overlap_nbases_min = 12
        String quant_mode = "TranscriptomeSAM GeneCounts"
        String read_files_command = "gunzip -c"
        String twopass_mode = "None"
        Int win_anchor_multimap_nmax = 50
        String extra_star_args = ""

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
        elif [[ "~{input_file_type}" == "FASTQ" ]]; then
            FASTQS_OPTION="~{sep=' ' fastqs}"
        fi

        echo "Extracting STAR index"
        mkdir -p "star_index"
        tar -C "star_index" -xzf "~{star_index}" --strip-components=1

        echo "Running STAR"
        STAR \
            --alignInsertionFlush "~{align_insertion_flush}" \
            --alignIntronMax ~{align_intron_max} \
            --alignMatesGapMax ~{align_mates_gap_max} \
            --alignSJDBoverhangMin ~{align_sjdb_overhang_min} \
            --alignSJstitchMismatchNmax "~{align_sj_stitch_mismatch_nmax}" \
            --alignSplicedMateMapLmin ~{align_spliced_mate_map_lmin} \
            --alignSplicedMateMapLminOverLmate ~{align_spliced_mate_map_lmin_over_lmate} \
            --chimJunctionOverhangMin ~{chim_junction_overhang_min} \
            --chimMultimapNmax ~{chim_multimap_nmax} \
            --chimMultimapScoreRange ~{chim_multimap_score_range} \
            --chimNonchimScoreDropMin ~{chim_nonchim_score_drop_min} \
            --chimOutJunctionFormat ~{chim_out_junction_format} \
            --chimOutType "~{chim_out_type}" \
            --chimScoreJunctionNonGTAG ~{chim_out_type} \
            --chimSegmentMin ~{chim_segment_min} \
            --genomeDir "~{genome_dir}" \
            --genomeLoad "~{genome_load}" \
            --limitSjdbInsertNsj ~{limit_sjdb_insert_nsj} \
            --outFilterIntronMotifs "~{out_filter_intron_motifs}" \
            --outFilterMatchNmin ~{out_filter_match_nmin} \
            --outFilterMatchNminOverLread ~{out_filter_match_nmin_over_lread} \
            --outFilterMismatchNmax ~{out_filter_mismatch_nmax} \
            --outFilterMismatchNoverLmax ~{out_filter_mismatch_nover_lmax} \
            --outFilterMultimapNmax ~{out_filter_multimap_nmax} \
            --outFilterScoreMinOverLread ~{out_filter_score_min_over_lread} \
            --outFilterType "~{out_filter_type}" \
            --outReadsUnmapped "~{out_reads_unmapped}" \
            --outSAMattrRGline "~{out_sam_attr_rg_line}" \
            --outSAMattributes "~{out_sam_attributes} \
            --outSAMstrandField "~{out_sam_strand_field}" \
            --outSAMtype "~{out_sam_type}" \
            --outSAMunmapped "~{out_sam_unmapped}" \
            --peOverlapMMp ~{pe_overlap_mmp} \
            --peOverlapNbasesMin ~{pe_overlap_nbases_min} \
            --quantMode "~{quant_mode}" \
            --readFilesCommand "~{read_files_command}" \
            --readFilesIn $FASTQS_OPTION \
            --runThreadN "~{n_threads} \
            --twopassMode "~{twopass_mode}" \
            --winAnchorMultimapNmax ~{win_anchor_multimap_nmax} \
            "~{extra_star_args}"

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
