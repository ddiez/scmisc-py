def sample_cells(adata, n, groupby=None, random_state=None):
  import scanpy as sc
  
  if groupby is None:
    return sc.pp.sample(adata, n=n, rng=random_state, copy=True)
  
  groups = adata.obs[groupby].unique()
  adata_list = []
  for group in groups:
    nsample = n
    tmp = adata[adata.obs[groupby] == group, :]
    if nsample > tmp.n_obs:
      nsample = tmp.n_obs
    adata_list.append(sc.pp.sample(tmp, n=nsample, rng=random_state, copy=True))
  
  return sc.concat(adata_list)