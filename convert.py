import pandas as pd
import numpy as np


def clean_up(vcf_in, txt_out, verbose=True):
    """Convert the whole genome VCF file to Genetic Lifehacks format and save to txt file."""

    # Reading VCF data and processing to get the Genetic Lifehacks format.
    df = pd.read_csv(
        vcf_in,
        compression="gzip",
        comment="#",
        sep="\t",
        names=[
            "chromosome",
            "position",
            "rsid",
            "reference",
            "alternative",
            "qual",
            "filter",
            "info",
            "format",
            "sample",
        ],
    )
    df = df.drop(columns=["qual", "filter", "format", "sample"])
    df["count"] = df["info"].str[3]
    df = df.drop(columns=["info"])
    df["genotype"] = np.where(
        df["count"] == "2",
        (df["alternative"] + df["alternative"]),
        (df["reference"] + df["alternative"]),
    )
    df = df.loc[
        :,
        [
            "rsid",
            "chromosome",
            "position",
            "genotype",
            "reference",
            "alternative",
            "count",
        ],
    ]

    """ Preview for quality checking if the the genotype column is constructed
    correctly. Rules: If count column == 2, the genotype is homozygous
    (genotype column should be double alternative column value),
    if count == 1 then heterozygous (genotype column should be reference + alternative columns) """
    if verbose:
        print(f"Preview top 40 rows:")
        print(df.head(40))

    # Dropping the remaining unnecessary columns and savint txt file.
    df = df.drop(columns=["reference", "alternative", "count"])
    df.to_csv(txt_out, sep="\t", mode="a", header=True, index=False)
    print(f"TXT file saved to {txt_out}")


if __name__ == "__main__":
    vcf_in = "annotated.snp.vcf.gz"
    txt_out = "glh_output.txt"
    clean_up(vcf_in, txt_out)
