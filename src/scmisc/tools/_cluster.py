def subcluster_group(adata, var, group, layer, resolution=0.4):
  import scanpy as sc
  cells = adata.obs[var] == group
  tmp = adata[cells, :].copy()
  sc.pp.pca(tmp)
  sc.pp.neighbors(tmp, n_pcs=10)
  sc.tl.leiden(tmp, key_added="_subcluster", resolution=resolution)

  adata.obs[var] = adata.obs[var].to_list()
  adata.obs[var][cells] = adata.obs[var][cells] + "_" + tmp.obs["_subcluster"].to_list()