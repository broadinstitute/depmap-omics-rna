version 1.0

import "../align_rna_sample/align_rna_sample.wdl" as align_rna_sample
import "../call_sr_rna_fusions/call_sr_rna_fusions.wdl" as call_sr_rna_fusions
import "../quantify_sr_rna/quantify_sr_rna.wdl" as quantify_sr_rna

workflow process_sr_rna {
    input {
        String sample_id
        String input_file_type
        Array[File]? fastqs
        File? cram_bam
        File? crai_bai
        File ref_fasta
        File? ref_fasta_index
        File star_index
        File gtf
        File targets
        Boolean stranded
        String? salmon_lib_type_override
    }

    if (input_file_type == "FASTQ") {
        call align_rna_sample.align_with_star as align_fastq_with_star {
            input:
                sample_id = sample_id,
                input_file_type = input_file_type,
                fastqs = fastqs,
                star_index = star_index
        }
    }

    if (input_file_type == "BAM") {
        call align_rna_sample.align_with_star as align_bam_with_star {
            input:
                sample_id = sample_id,
                input_file_type = input_file_type,
                cram_bam = cram_bam,
                star_index = star_index
        }
    }

    if (input_file_type == "CRAM") {
        call align_rna_sample.align_with_star as align_cram_with_star {
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

    call call_sr_rna_fusions.call_fusions_with_arriba as call_fusions_with_arriba {
        input:
            sample_id = sample_id,
            analysis_ready_bam = select_first([
                align_bam_with_star.analysis_ready_bam,
                align_cram_with_star.analysis_ready_bam
            ]),
            ref_fasta = ref_fasta,
            gtf = gtf
    }

    call quantify_sr_rna.quantify_with_salmon as quantify_with_salmon {
        input:
            sample_id = sample_id,
            transcriptome_bam = select_first([
                align_bam_with_star.transcriptome_bam,
                align_cram_with_star.transcriptome_bam
            ]),
            gtf = gtf,
            targets = targets,
            stranded = stranded,
            lib_type_override = salmon_lib_type_override
    }

    output {
        File reads_per_gene = select_first([
            align_fastq_with_star.reads_per_gene,
            align_bam_with_star.reads_per_gene,
            align_cram_with_star.reads_per_gene
        ])
        File junctions = select_first([
            align_fastq_with_star.junctions,
            align_bam_with_star.junctions,
            align_cram_with_star.junctions
        ])
        File fusions = call_fusions_with_arriba.fusions
        File fusions_discarded = call_fusions_with_arriba.fusions_discarded
        File quant_transcripts = quantify_with_salmon.quant_transcripts
        File quant_genes = quantify_with_salmon.quant_genes
    }
}
