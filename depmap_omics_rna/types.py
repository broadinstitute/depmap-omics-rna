from typing import Optional, TypeVar

import httpx
import pandas as pd
import pandera as pa
import pandera.typing
from nebelung.types import CoercedDataFrame
from pandera.typing import Series
from pydantic import BaseModel

from gumbo_gql_client.gumbo_client import GumboClient as AriadneGumboClient


class GumboClient(AriadneGumboClient):
    def __init__(
        self,
        url: str,
        username: str,
        headers: dict[str, str],
        http_client: Optional[httpx.Client] = None,
    ):
        super().__init__(url=url, headers=headers, http_client=http_client)
        self.username = username  # store username on this object for use in mutations


class GumboAlignment(CoercedDataFrame):
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    omics_sequencing_id: Series[pd.StringDtype]
    sequencing_alignment_source: Series[pd.StringDtype] = pa.Field(isin=["GP", "CDS"])
    reference_genome: Series[pd.StringDtype]
    url: Series[pd.StringDtype] = pa.Field(unique=True)
    index_url: Series[pd.StringDtype] = pa.Field(nullable=True)
    size: Series[pd.Int64Dtype] = pa.Field(unique=True)
    stranded: Series[pd.BooleanDtype]


class TerraSample(CoercedDataFrame):
    sample_id: Series[pd.StringDtype] = pa.Field(unique=True)
    model_id: Series[pd.StringDtype]
    cell_line_name: Series[pd.StringDtype]
    stripped_cell_line_name: Series[pd.StringDtype]
    model_condition_id: Series[pd.StringDtype]
    omics_profile_id: Series[pd.StringDtype]
    delivery_cram_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    delivery_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    delivery_file_format: Series[pd.StringDtype] = pa.Field(
        isin={"CRAM", "BAM"}, nullable=True
    )
    delivery_cram_bam_size: Series[pd.Int64Dtype] = pa.Field(unique=True)
    delivery_stranded: Series[pd.BooleanDtype]
    delivery_ref: Series[pd.StringDtype]
    delivery_ref_fasta: Series[pd.StringDtype]
    delivery_ref_fasta_index: Series[pd.StringDtype]


class GcsObject(CoercedDataFrame):
    url: Series[pd.StringDtype]
    size: Series[pd.Int64Dtype]
    crc32c_hash: Series[pd.StringDtype]
    created_at: Series[pd.Timestamp]


class DeltaJob(BaseModel):
    workflow_name: str
    entity_type: str
    entity_set_type: str
    entity_id_col: str
    expression: str
    input_cols: set[str] | None = None
    output_cols: set[str] | None = None
    resubmit_n_times: int = 0
    force_retry: bool = False
    use_callcache: bool = True
    use_reference_disks: bool = False
    memory_retry_multiplier: float = 1.0
    max_n_entities: int | None = None
    dry_run: bool = False


PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = pandera.typing.DataFrame
T = TypeVar("T")
