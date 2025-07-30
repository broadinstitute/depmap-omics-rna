import pandas as pd
from nebelung.terra_workspace import TerraWorkspace

from depmap_omics_rna.utils import get_gcs_object_metadata

pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.max_info_columns", 30)
pd.set_option("display.max_info_rows", 20)
pd.set_option("display.max_rows", 20)
pd.set_option("display.max_seq_items", None)
pd.set_option("display.width", 200)
pd.set_option("expand_frame_repr", True)
pd.set_option("mode.chained_assignment", "warn")

ws = TerraWorkspace("broad-firecloud-ccle", "depmap-omics-rna")
samples = ws.get_entities("sample")

df = (
    samples.loc[
        :,
        ["sample_id", "reads_per_gene", "star_junctions"],
    ]
    .melt(id_vars="sample_id", var_name="col", value_name="url")
    .dropna()
)

blobs = get_gcs_object_metadata(df["url"], "depmap-omics")

df = df.merge(blobs, on="url", how="left")

df = df.sort_values("size")

ws.upload_entities(samples.loc[samples["sample_id"].isin(["CDS-sFk4sv", "CDS-Roe0EA"])])
