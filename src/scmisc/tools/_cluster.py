def subcluster_group(adata, reduction_method, source_key, group, target_key, layer=None, **kargv):
  import numpy as np

  # If the target key exists, convert it to string to allow concatenation later.
  if np.any(adata.obs.columns == target_key):
    adata.obs[target_key] = adata.obs[target_key].astype(str)

  # If the target key does not exist, fill it with the source key but convert to strings.
  if ~np.any(adata.obs.columns == target_key):
    adata.obs[target_key] = adata.obs[source_key].astype(str)

  if not isinstance(group, list):
     group = list(group)
  
  cells = adata.obs_names[adata.obs[source_key].isin(group)].to_list()
  tmp = adata[cells, :].copy()

  if layer is not None:
    tmp.X = tmp.layers[layer]
  
  match reduction_method:
    case "pca":
        subclusters = subcluster_group_pca(tmp, **kargv)
    
    case "harmony":
        subclusters = subcluster_group_harmony(tmp, **kargv)
    
    case "scvi":
        subclusters = subcluster_group_scvi(tmp, **kargv)
  
  adata.obs.loc[cells, target_key] = adata.obs.loc[cells, target_key] + "_" + subclusters.astype(str)

def subcluster_group_pca(adata, n_pcs=10, resolution=0.4):
  import scanpy as sc
  import numpy as np

  sc.pp.pca(adata)
  
  sc.pp.neighbors(adata, n_pcs=n_pcs)
  sc.tl.leiden(adata, key_added=".subcluster", resolution=resolution, flavor="igraph", n_iterations=2, directed=False)

  return adata.obs[".subcluster"]

def subcluster_group_harmony(adata, resolution=0.4):
  # TODO: add support
  import warnings
  warnings.warn("Not supported yet")

def subcluster_group_scvi(adata, layer=None, batch_key=None, categorical_covariate_keys=None, seed=0, resolution=0.4):
  import scanpy as sc
  import numpy as np
  import scvi

  scvi.settings.seed = seed

  scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key, categorical_covariate_keys=categorical_covariate_keys)
  model = scvi.model.SCVI(adata)
  model.train(max_epochs=400, early_stopping=True, early_stopping_monitor="validation_loss")
  adata.obsm["X_scvi"] = model.get_latent_representation()

  sc.pp.neighbors(adata, use_rep="X_scvi")
  sc.tl.leiden(adata, key_added=".subcluster", resolution=resolution, flavor="igraph", n_iterations=2, directed=False)

  return adata.obs[".subcluster"]