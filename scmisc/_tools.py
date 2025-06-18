def subcluster_group_pca(adata, source_key, group, target_key, layer=None, n_pcs=10, resolution=0.4, seed=0):
  import scanpy as sc
  import numpy as np

  target_key_tmp = "." + target_key

  cells = adata.obs_names[adata.obs[source_key] == group].to_list()
  tmp = adata[cells, :].copy()

  if layer is not None:
    tmp.X = tmp.layers[layer]

  sc.pp.pca(tmp)
  
  sc.pp.neighbors(tmp, n_pcs=n_pcs)
  sc.tl.leiden(tmp, key_added=target_key_tmp, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)
  
  if np.any(adata.obs.columns == target_key):
     adata.obs[target_key] = adata.obs[target_key].to_list()
  if ~np.any(adata.obs.columns == target_key):
      adata.obs[target_key] = adata.obs[source_key].to_list()
    
  adata.obs.loc[cells, target_key] = adata.obs.loc[cells, target_key] + "_" + tmp.obs[target_key_tmp].to_list()

def subcluster_group_scvi(adata, source_key, group, target_key, layer="counts", resolution=0.4, batch_key=None, categorical_covariate_keys=None, seed=0):
  import scanpy as sc
  import numpy as np
  import scvi

  scvi.settings.seed = seed

  target_key_tmp = "." + target_key

  cells = adata.obs_names[adata.obs[source_key] == group].to_list()
  tmp = adata[cells, :].copy()

  scvi.model.SCVI.setup_anndata(tmp, layer=layer, batch_key=batch_key, categorical_covariate_keys=categorical_covariate_keys)
  model = scvi.model.SCVI(tmp)
  model.train(max_epochs=400, early_stopping=True, early_stopping_monitor="validation_loss")
  tmp.obsm["X_scvi"] = model.get_latent_representation()

  sc.pp.neighbors(tmp, use_rep="X_scvi")
  sc.tl.leiden(tmp, key_added=target_key_tmp, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)

  if np.any(adata.obs.columns == target_key):
     adata.obs[target_key] = adata.obs[target_key].to_list()
  if ~np.any(adata.obs.columns == target_key):
      adata.obs[target_key] = adata.obs[source_key].to_list()
    
  adata.obs.loc[cells, target_key] = adata.obs.loc[cells, target_key] + "_" + tmp.obs[target_key_tmp].to_list()

def filter_cells_by_genes(adata, genes, cutoff=0, layer="counts"):
  for gene in genes:
    if gene in adata.var_names:
      sel = adata[:, adata.var_names == gene].layers[layer].toarray() > cutoff
      adata = adata[~sel, :]
  return adata

def calculate_dims(ngroups, nrows=None, ncols=None):
  import numpy as np

  if ncols is None and nrows is None:
      ncols = np.int32(np.ceil(np.sqrt(ngroups)))
      nrows = np.int32(np.ceil(ngroups/ncols))

  if ncols is None and nrows is not None:
    ncols = np.int32(np.ceil(ngroups/nrows))

  if nrows is None and ncols is not None:
    nrows = np.int32(np.ceil(ngroups/ncols))

  return nrows, ncols
