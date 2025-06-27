from typing import Iterable, Union

import pandas as pd
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


def anti_join(
    x: pd.DataFrame, y: pd.DataFrame, on: Union[str, Iterable[str]]
) -> pd.DataFrame:
    """
    Anti join two data frames.

    :param x: a base data frame
    :param y: a data frame to use for filtering
    :param on: the columns to anti-join on
    :return: a data frame
    """

    if len(y) == 0:
        return x

    # make a data frame of just the join columns
    dummy = y.loc[:, on]  # pyright: ignore

    # convert to data frame if `on` was a single column
    if isinstance(dummy, pd.Series):
        dummy = dummy.to_frame()

    dummy.loc[:, "dummy_col"] = 1  # indicator variable

    # attempt to join left data frame (`x`) to the dummy data frame
    merged = x.merge(dummy, on=on, how="left")

    # keep only the non-matches
    return merged.loc[merged["dummy_col"].isna(), x.columns.tolist()]


ws = TerraWorkspace("broad-firecloud-ccle", "depmap-omics-rna-dev")
samples = ws.get_entities("sample")

samples.dropna(subset=["star_junctions_new"], inplace=True)

olds = []

for _, r in samples.iterrows():
    olds.append(
        pd.read_csv(
            r["star_junctions"],
            sep="\t",
            names=[
                "chr",
                "intron_start",
                "intron_end",
                "strand",
                "motif",
                "annotated",
                "unique_reads",
                "multi_reads",
                "max_overhang",
                "unique_read_id",
                "donor_overhang",
                "acceptor_overhang",
            ],
        ).assign(sample_id=r["sample_id"])
    )

old = pd.concat(olds)
old.to_parquet("./data/old.parquet", index=False)

news = []

for _, r in samples.iterrows():
    news.append(
        pd.read_csv(
            r["star_junctions_new"],
            sep="\t",
            names=[
                "chr",
                "intron_start",
                "intron_end",
                "strand",
                "motif",
                "annotated",
                "unique_reads",
                "multi_reads",
                "max_overhang",
                "unique_read_id",
                "donor_overhang",
                "acceptor_overhang",
            ],
        ).assign(sample_id=r["sample_id"])
    )

new = pd.concat(news)
new.to_parquet("./data/new.parquet", index=False)

anti_join(
    old,
    new,
    on=[
        "chr",
        "intron_start",
        "intron_end",
        "strand",
        "motif",
        "annotated",
        "unique_reads",
        "multi_reads",
        "max_overhang",
        "unique_read_id",
        "donor_overhang",
        "acceptor_overhang",
    ],
)

anti_join(old, new, on=["chr", "intron_start", "intron_end"])

anti_join(
    new,
    old,
    on=[
        "chr",
        "intron_start",
        "intron_end",
        "strand",
        "motif",
        "annotated",
        "unique_reads",
        "multi_reads",
        "max_overhang",
        "unique_read_id",
        "donor_overhang",
        "acceptor_overhang",
    ],
)

anti_join(new, old, on=["chr", "intron_start", "intron_end"])
