from functools import partial

import pandas as pd
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import type_data_frame
from pd_flatten import pd_flatten

from depmap_omics_rna.types import GumboAlignment, GumboClient, TerraSample
from depmap_omics_rna.utils import model_to_df


def refresh_terra_samples(
    terra_workspace: TerraWorkspace,
    gumbo_client: GumboClient,
    ref_urls: dict[str, dict[str, str]],
) -> None:
    """
    Update the Terra `sample` data table using ground truth data from Gumbo.

    :param terra_workspace: a `TerraWorkspace` instance
    :param gumbo_client: a `GumboClient` instance
    :param ref_urls: a nested dictionary of genomes and their reference file URLs
    """

    # get long data frame of both GP-delivered and CDS (analysis ready) CRAM/BAMs
    alignments = model_to_df(
        gumbo_client.rna_sequencing_alignments(),
        GumboAlignment,
        mutator=partial(pd_flatten, name_columns_with_parent=False),
    )

    # make wide, separating delivery and old analysis-ready CRAM/BAMs
    samples = (
        alignments.loc[alignments["sequencing_alignment_source"].eq("GP")]
        .drop(columns=["sequencing_alignment_source"])
        .rename(
            columns={
                "omics_sequencing_id": "sample_id",
                "sequencing_alignment_id": "delivery_sequencing_alignment_id",
                "url": "delivery_cram_bam",
                "size": "delivery_cram_bam_size",
                "index_url": "delivery_crai_bai",
                "reference_genome": "delivery_ref",
                "stranded": "delivery_stranded",
            }
        )
        .merge(
            alignments.loc[
                alignments["sequencing_alignment_source"].eq("CDS")
                & alignments["reference_genome"].eq("hg38")
            ]
            .drop(columns=["sequencing_alignment_source"])
            .rename(
                columns={
                    "omics_sequencing_id": "sample_id",
                    "sequencing_alignment_id": "old_analysis_ready_sequencing_alignment_id",
                    "url": "old_analysis_ready_bam",
                    "size": "old_analysis_ready_bam_size",
                    "index_url": "old_analysis_ready_bai",
                    "reference_genome": "old_analysis_ready_ref",
                    "stranded": "old_analysis_ready_stranded",
                }
            ),
            how="outer",
            on=[
                "sample_id",
                "model_id",
                "model_condition_id",
                "omics_profile_id",
                "cell_line_name",
                "stripped_cell_line_name",
            ],
        )
    )

    # use old analysis-ready BAMs as a backup for the original GP-delivered ones
    samples["delivery_cram_bam"] = samples["delivery_cram_bam"].fillna(
        samples["old_analysis_ready_bam"]
    )

    samples["delivery_cram_bam_size"] = samples["delivery_cram_bam_size"].fillna(
        samples["old_analysis_ready_bam_size"]
    )

    samples["delivery_crai_bai"] = samples["delivery_crai_bai"].fillna(
        samples["old_analysis_ready_bai"]
    )

    samples["delivery_ref"] = samples["delivery_ref"].fillna(
        samples["old_analysis_ready_ref"]
    )

    samples["delivery_stranded"] = samples["delivery_stranded"].fillna(
        samples["old_analysis_ready_stranded"]
    )

    samples["delivery_file_format"] = (
        samples["delivery_cram_bam"].str.rsplit(".", n=1).str.get(1).str.upper()
    )

    # set reference genome columns
    samples = set_ref_urls(samples, ref_urls)

    # validate types
    samples = type_data_frame(samples, TerraSample, remove_unknown_cols=True)

    # delete obsolete samples (e.g. ones that have been blacklisted since the last sync)
    terra_samples = terra_workspace.get_entities("sample")

    if len(terra_samples) > 0:
        terra_workspace.delete_entities(
            entity_type="sample",
            entity_ids=set(terra_samples["sample_id"]).difference(
                set(samples["sample_id"])
            ),
        )

    sample_ids = samples.pop("sample_id")
    samples.insert(0, "entity:sample_id", sample_ids)
    terra_workspace.upload_entities(df=samples)


def set_ref_urls(
    samples: pd.DataFrame, ref_urls: dict[str, dict[str, str]]
) -> pd.DataFrame:
    """
    Populate columns in the `samples` data frame with URLs for reference genome files.

    :param samples: the data frame of sample data
    :param ref_urls: a dictionary of reference genome names (e.g. "hg38") and URLs of
    files
    :return: the `samples` data frame with columns for reference genome URLs
    """

    # join reference URLs to for `delivery_` and `analysis_ready_` CRAM/BAMs
    ref_df = pd.DataFrame(ref_urls.values())
    ref_df["ref"] = ref_urls.keys()
    ref_df.columns = "delivery_" + ref_df.columns

    return samples.merge(ref_df, how="left", on="delivery_ref")
