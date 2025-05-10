import ast
import pandas as pd
from datasets import load_dataset

# read in edist matrix from tahoe perturbations
edist = pd.read_table("Tahoe96M_edist_matrix_within_plate.tsv", sep="\t")

# parse out the drug names from the conceentrations in the column names that aren't cell_line
# make a new df with 1 column per cell line, drug_name, drug_concentration, e_distance (which are the values in the matrix)
edist_long = edist.melt(id_vars=["cell_line"], var_name="drug_info", value_name="e_distance")
# make the drug_info column into a tuple
edist_long["drug_info"] = edist_long["drug_info"].apply(lambda x: ast.literal_eval(x)[0])
# make new columns for drug_name, concentration, unit
edist_long[["drug_name", "concentration", "unit"]] = pd.DataFrame(
    edist_long["drug_info"].tolist(), index=edist_long.index
)


# load tahoe drug metadata to add on targets, MoA, clinical trial info, pubchem id, etc.
ds = load_dataset("tahoebio/Tahoe-100M", "drug_metadata")
drug_metadata = ds["train"].to_pandas()
drug_name = drug_metadata["drug"]

# load tahoe cell line metadata to get cell line names, cell line ids, and cell line types
cs = load_dataset("tahoebio/Tahoe-100M", "cell_line_metadata")
cell_line_metadata = cs["train"].to_pandas()
cell_line = cell_line_metadata["Cell_ID_Cellosaur"]
# optionally listify all the unique values in columns in cell_line_metadata that aren't cell_name e.g. such that Driver_Gene_Symbol looks like ["CDKN2A", "CDKN2B", "KRAS","SMARCA4","STK11"] for A549 cell_name

# merge the drug and cell line metadata with the edist_long df
edist_long_drug_metadata = edist_long.merge(drug_metadata, left_on="drug_name", right_on="drug", how="inner")
edist_long_drug_metadata.shape

# for each unique drug_name and cell_line, select the concentration that has the largest e_distance
top_edist_per_drug_cell_line_combo = (
    edist_long_drug_metadata.sort_values("e_distance", ascending=False)
    .drop_duplicates(subset=["drug_name", "cell_line"])
    .reset_index(drop=True)
)
top_edist_per_drug_cell_line_combo.to_csv("tahoe_largest_edist_per_drug_cell_line_combo.csv", index=False)
