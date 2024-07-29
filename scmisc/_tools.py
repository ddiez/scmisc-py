def subcluster_group(adata, var, group, layer, resolution=0.4):
  import scanpy as sc
  cells = adata.obs[var] == group
  tmp = adata[cells, :].copy()
  sc.pp.pca(tmp)
  sc.pp.neighbors(tmp, n_pcs=10)
  sc.tl.leiden(tmp, key_added="_subcluster", resolution=resolution)

  adata.obs[var] = adata.obs[var].to_list()
  adata.obs[var][cells] = adata.obs[var][cells] + "_" + tmp.obs["_subcluster"].to_list()

def calculate_dims(ngroups, nrows=None, ncols=None):
  import numpy as np

  if ncols is None and nrows is None:
      ncols = np.int32(np.ceil(np.sqrt(ngroups)))
      nrows = int(np.ceil(ngroups/ncols))

  if ncols is None and nrows is not None:
    ncols = int(np.ceil(ngroups/nrows))

  if nrows is None and ncols is not None:
    nrows = int(np.ceil(ngroups/ncols))

  return nrows, ncols
