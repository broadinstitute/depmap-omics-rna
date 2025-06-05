import firecloud_api_cds.api as firecloud_api
import pandas as pd
from google.cloud import storage
from nebelung.terra_workspace import TerraWorkspace
from nebelung.utils import batch_evenly, call_firecloud_api

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

old_ws = TerraWorkspace("broad-firecloud-ccle", "DepMap_hg38_RNAseq")
old_samples = (
    old_ws.get_entities("sample").loc[:, ["sample_id", "star_junctions"]].dropna()
)

ws = TerraWorkspace("broad-firecloud-ccle", "depmap-omics-rna")
samples = ws.get_entities("sample")

old_samples = old_samples.loc[old_samples["sample_id"].isin(samples["sample_id"])]

src_bucket_name = call_firecloud_api(
    firecloud_api.get_workspace,
    namespace=old_ws.workspace_namespace,
    workspace=old_ws.workspace_name,
    fields=["workspace.bucketName"],
)["workspace"]["bucketName"]

dest_bucket_name = call_firecloud_api(
    firecloud_api.get_workspace,
    namespace=ws.workspace_namespace,
    workspace=ws.workspace_name,
    fields=["workspace.bucketName"],
)["workspace"]["bucketName"]

storage_client = storage.Client()
src_bucket = storage_client.bucket(src_bucket_name)
dest_bucket = storage_client.bucket(dest_bucket_name)

urls = old_samples.copy().rename(columns={"star_junctions": "url"})
urls["new_url"] = (
    "gs://" + dest_bucket.name + "/legacy/" + urls["sample_id"] + ".SJ.out.tab.gz"
)

for batch in batch_evenly(urls.sort_index(ascending=False), 500):
    with storage_client.batch(raise_exception=False):
        for _, r in batch.iterrows():
            source_blob = storage.Blob.from_string(r["url"], client=storage_client)
            destination_blob = storage.Blob.from_string(
                r["new_url"], client=storage_client
            )

            blob_copy = src_bucket.copy_blob(
                source_blob, dest_bucket, destination_blob.name
            )
            print(blob_copy.name)

ws.upload_entities(
    urls.rename(columns={"new_url": "star_junctions"}).drop(columns=["url"])
)
