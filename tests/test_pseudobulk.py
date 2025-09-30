import scmisc
import scanpy as sc
import numpy as np
import pandas as pd


def test_pseudobulk():
  X = np.arange(1, 17).reshape(-1, 2)
  adata = sc.AnnData(X)
  adata.obs_names = [ "cell_" + str(i+1) for i in range(8) ]
  adata.var_names = [ "gene_" + str(i+1) for i in range(2) ]
  adata.obs["sample"] = np.repeat(["D1", "D2"], repeats=4)
  adata.obs["cluster"] = ["C1", "C2", "C3", "C1", "C2", "C1", "C1", "C2"]
  meta = pd.DataFrame({
    "Group1": ["F", "F", "M", "M"],
    "Group2": ["H", "D", "H", "D"]},
    index=["D1", "D2", "D3", "D4"]
  )

  mdata = scmisc.tools.pseudobulk(adata, groupby="sample", splitby="cluster", aggregate="sum", metadata=None, return_mudata=False)

  assert mdata
  assert mdata["C1"] is not None
  assert mdata["C2"] is not None
  assert mdata["C3"] is not None

  print(mdata["C1"])
  print(mdata["C2"])
  print(mdata["C3"])

  assert (mdata["C1"].loc["D1", "gene_1"] == 8)
  assert (mdata["C1"].loc["D1", "gene_2"] == 10)
  assert (mdata["C1"].loc["D2", "gene_1"] == 24)
  assert (mdata["C1"].loc["D2", "gene_2"] == 26)

  assert (mdata["C2"].loc["D1", "gene_1"] == 3)
  assert (mdata["C2"].loc["D1", "gene_2"] == 4)
  assert (mdata["C2"].loc["D2", "gene_1"] == 24)
  assert (mdata["C2"].loc["D2", "gene_2"] == 26)

  assert (mdata["C3"].loc["D1", "gene_1"] == 5)
  assert (mdata["C3"].loc["D1", "gene_2"] == 6)
  assert (mdata["C3"].loc["D2", "gene_1"] == 0)
  assert (mdata["C3"].loc["D2", "gene_2"] == 0)

