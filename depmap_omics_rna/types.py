from typing import Any, Optional, TypeVar

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
    omics_sequencing_id: Series[pd.StringDtype]
    sequencing_alignment_source: Series[pd.StringDtype] = pa.Field(isin=["GP", "CDS"])
    url: Series[pd.StringDtype] = pa.Field(unique=True)
    index_url: Series[pd.StringDtype] = pa.Field(unique=True)


class TerraSample(CoercedDataFrame):
    sample_id: Series[pd.StringDtype] = pa.Field(unique=True)
    delivery_cram_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    delivery_crai_bai: Series[pd.StringDtype] = pa.Field(nullable=True)
    delivery_file_format: Series[pd.StringDtype] = pa.Field(
        isin={"CRAM", "BAM"}, nullable=True
    )
    analysis_ready_bam: Series[pd.StringDtype] = pa.Field(nullable=True)
    analysis_ready_bai: Series[pd.StringDtype] = pa.Field(nullable=True)


class GumboTaskEntity(CoercedDataFrame):
    id: Series[pd.StringDtype]
    sequencing_id: Series[pd.StringDtype] = pa.Field(nullable=True)


class GumboTaskResult(CoercedDataFrame):
    sample_id: Series[pd.StringDtype]
    label: Series[pd.StringDtype]
    url: Series[pd.StringDtype] = pa.Field(nullable=True)
    value: Series[dict[str, Any]] = pa.Field(nullable=True)
    workflow_name: Series[pd.StringDtype] = pa.Field(nullable=True)
    workflow_version: Series[pd.StringDtype] = pa.Field(nullable=True)
    completed_at: Series[pd.Timestamp]


class GcsObject(CoercedDataFrame):
    url: Series[pd.StringDtype]
    size: Series[pd.Int64Dtype]
    crc32c_hash: Series[pd.StringDtype]
    created_at: Series[pd.Timestamp]


PydanticBaseModel = TypeVar("PydanticBaseModel", bound=BaseModel)
PanderaBaseSchema = TypeVar("PanderaBaseSchema", bound=CoercedDataFrame)
TypedDataFrame = pandera.typing.DataFrame
T = TypeVar("T")
