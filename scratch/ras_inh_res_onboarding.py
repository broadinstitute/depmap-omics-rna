import base64
import hashlib
import json

import pandas as pd
from google.cloud import storage
from google.cloud.storage import Blob
from nebelung.terra_workspace import TerraWorkspace

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

storage_client = storage.Client()


def list_blobs(
    bucket_name: str, prefix: str | None = None, glob: str | None = None
) -> pd.DataFrame:
    """
    Get the names and sizes of existing blobs in a GCS bucket.

    :param bucket_name: the name of the GCS bucket
    :param prefix: an optional prefix for listing
    :param glob: an optional glob for listing
    :return: a data frame of object names and sizes
    """

    if prefix is not None and glob is not None:
        raise ValueError("At most one of `prefix` and `glob` can be specified")
    elif prefix is not None:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name, prefix=prefix, delimiter="/"
        ).pages
    elif glob is not None:
        pages = storage_client.list_blobs(
            bucket_or_name=bucket_name, match_glob=glob
        ).pages
    else:
        pages = storage_client.list_blobs(bucket_or_name=bucket_name).pages

    blobs = []

    for page in pages:
        blobs.extend(
            [
                {
                    "url": "gs://" + bucket_name + "/" + x.name,
                    "md5": base64.b64decode(x.md5_hash).hex(),
                    "size": x.size,
                    "gcs_obj_updated_at": x.updated,
                }
                for x in page
            ]
        )

    return pd.DataFrame(blobs)


def get_hash(url: str) -> str:
    print(url)
    return Blob.from_string(url).download_as_text(client=storage_client).split(" ")[0]


blobs = list_blobs(bucket_name="depmap-moffitt", glob="*.fastq.*")
blobs["sample_id"] = blobs["url"].str.extract(r"^.+/(.+)_\d\.fastq\..+$")

fastqs = (
    blobs.loc[blobs["url"].str.endswith(".gz")]
    .copy()
    .rename(columns={"url": "fastq_url"})
)

hashes = (
    blobs.loc[blobs["url"].str.endswith(".md5")]
    .copy()
    .rename(columns={"url": "md5_url"})
)

hashes["fastq_url"] = hashes["md5_url"].str.removesuffix(".md5")

fastqs = fastqs[["sample_id", "fastq_url", "md5"]].merge(
    hashes[["fastq_url", "md5_url"]],
    on="fastq_url",
    how="inner",
    validate="one_to_one",
)

fastqs["md5_check"] = fastqs["md5_url"].apply(get_hash)

assert bool(fastqs["md5"].eq(fastqs["md5_check"]).all())
assert bool(~fastqs["md5"].duplicated().any())

samples = (
    fastqs.sort_values("fastq_url")
    .groupby("sample_id")["fastq_url"]
    .agg(list)
    .reset_index()
    .rename(columns={"fastq_url": "fastqs"})
)

samples["fastqs"] = samples["fastqs"].apply(json.dumps)

metadata = pd.read_csv("../data/ras_inh_metadata.csv", dtype="string")
replicates = pd.read_csv("../data/ras_inh_replicates.csv").convert_dtypes()

samples = samples.merge(
    metadata, on="sample_id", how="inner", validate="one_to_one"
).merge(replicates, on="sample_id", how="inner", validate="one_to_one")

tw = TerraWorkspace("broad-firecloud-ccle", "depmap-ras-inhibitor-resistant")
tw.upload_entities(samples)
