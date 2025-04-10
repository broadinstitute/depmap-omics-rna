version 1.0

import "../align_rna_sample/align_rna_sample.wdl" as align_rna_sample
import "../call_sr_rna_fusions/call_sr_rna_fusions.wdl" as call_sr_rna_fusions
import "../quantify_sr_rna/quantify_sr_rna.wdl" as quantify_sr_rna

workflow process_sr_rna {
    input {
        String workflow_version = "1.0" # internal semver
        String workflow_source_url # populated automatically with URL of this script

        String sample_id
        String cram_or_bam
        File cram_bam
        File? crai_bai
        File ref_fasta
        File? ref_fasta_index
        File star_index
        File gtf
        File targets
        Boolean stranded
        String? salmon_lib_type_override
    }

    if (cram_or_bam == "BAM") {
        call align_rna_sample.align_with_star as align_bam_with_star {
            input:
                sample_id = sample_id,
                cram_or_bam = cram_or_bam,
                cram_bam = cram_bam,
                star_index = star_index
        }
    }

    if (cram_or_bam == "CRAM") {
        call align_rna_sample.align_with_star as align_cram_with_star {
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
            align_bam_with_star.reads_per_gene,
            align_cram_with_star.reads_per_gene
        ])
        File fusions = call_fusions_with_arriba.fusions
        File fusions_discarded = call_fusions_with_arriba.fusions_discarded
        File quant_transcripts = quantify_with_salmon.quant_transcripts
        File quant_genes = quantify_with_salmon.quant_genes
    }
}
