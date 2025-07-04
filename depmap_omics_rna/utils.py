import logging
from pathlib import Path
from typing import Any, Callable, Iterable, Type

import pandas as pd
from google.cloud import secretmanager_v1, storage
from nebelung.terra_workflow import TerraWorkflow
from nebelung.utils import batch_evenly, type_data_frame

from depmap_omics_rna.types import (
    GcsObject,
    PanderaBaseSchema,
    PydanticBaseModel,
    TypedDataFrame,
)


def get_hasura_creds(gumbo_env: str) -> dict[str, str]:
    """
    Get URL and password for Hasura GraphQL API.

    :param gumbo_env: the Gumbo env to get credentials for ('staging' or 'prod')
    :return: a dictionary with the GraphQL API URL and password
    """

    logging.info(f"Getting Hasura credentials for {gumbo_env}")

    return {
        "url": get_secret_from_sm(
            f"projects/814840278102/secrets/hasura-{gumbo_env}-api-url/versions/latest"
        ),
        "password": get_secret_from_sm(
            f"projects/814840278102/secrets/hasura-admin-secret-{gumbo_env}/versions/latest"
        ),
    }


def get_secret_from_sm(name: str) -> str:
    """
    Get the value of a secret from GCP Secret Manager.

    :param name: the fully-qualified name of a secret
    :return: the secret's decoded value
    """

    client = secretmanager_v1.SecretManagerServiceClient()
    request = secretmanager_v1.AccessSecretVersionRequest(mapping={"name": name})
    response = client.access_secret_version(request=request)
    return response.payload.data.decode()


def model_to_df(
    model: PydanticBaseModel,
    pandera_schema: Type[PanderaBaseSchema],
    records_key: str = "records",
    remove_unknown_cols: bool = False,
    mutator: Callable[[pd.DataFrame], pd.DataFrame] = lambda _: _,
) -> TypedDataFrame[PanderaBaseSchema]:
    """
    Dump a Pydantic model and convert it to a data frame typed by a Pandera schema.

    :param model: a Pydandict model containing a list of objects keyed by `records_key`
    :param pandera_schema: the Pandera schema to cast the model to
    :param records_key: the key/method name in `model` containing the records
    :param remove_unknown_cols: remove columns not specified in the schema
    :param mutator: an optional function to call on the data frame before typing (e.g.
    to rename columns to expected Pydantic field names)
    """

    records = model.model_dump()[records_key]

    df = pd.DataFrame(records)
    df = mutator(df)
    return type_data_frame(df, pandera_schema, remove_unknown_cols)


def get_gcs_object_metadata(
    urls: Iterable[str], gcp_project_id: str
) -> TypedDataFrame[GcsObject]:
    """
    Check existence and get metadata (size, hash, etc.) of GCS objects.

    :param urls: iterable of GCS URLs
    :param gcp_project_id: the ID of a GCP project to use for billing
    :return: data frame of object URLs and metadata
    """

    logging.info(f"Getting metadata about {len(list(urls))} GCS objects")
    storage_client = storage.Client(project=gcp_project_id)
    blobs = {}

    # the GCS batch context below has a max batch size of 1000, so do this outer layer
    # of batching, too)
    for batch in batch_evenly(urls, max_batch_size=200):
        with storage_client.batch(raise_exception=False):
            for url in batch:
                blob = storage.Blob.from_string(url, client=storage_client)
                bucket = storage_client.bucket(
                    blob.bucket.name, user_project=gcp_project_id
                )
                blob = bucket.get_blob(blob.name)
                blobs[url] = blob

    df = pd.DataFrame(
        [
            {
                "url": k,
                "size": v.size,
                "crc32c_hash": v.crc32c,
                "created_at": v.time_created,
            }
            for k, v in blobs.items()
        ]
    )

    # batching without raising exceptions makes all columns NA if an object was missing
    df = df.dropna()

    return type_data_frame(df, GcsObject)


def make_workflow_from_config(
    config: dict[str, Any], workflow_name: str, **kwargs: Any
) -> TerraWorkflow:
    """
    Make a TerraWorkflow object from a config entry.

    :param config: a config dictionary
    :param workflow_name: the name of the workflow referenced in the config
    :return: a TerraWorkflow instance
    """

    return TerraWorkflow(
        method_namespace=config["terra"][workflow_name]["method_namespace"],
        method_name=config["terra"][workflow_name]["method_name"],
        method_config_namespace=config["terra"][workflow_name][
            "method_config_namespace"
        ],
        method_config_name=config["terra"][workflow_name]["method_config_name"],
        method_synopsis=config["terra"][workflow_name]["method_synopsis"],
        workflow_wdl_path=Path(
            config["terra"][workflow_name]["workflow_wdl_path"]
        ).resolve(),
        method_config_json_path=Path(
            config["terra"][workflow_name]["method_config_json_path"]
        ).resolve(),
        **kwargs,
    )
