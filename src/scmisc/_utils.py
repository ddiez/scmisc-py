from anndata import AnnData
from mudata import MuData

def fetch_data(x, mod=None, obsm=None, obs=None, features=None, layer=None, use_raw=None):
    
    if isinstance(x, AnnData) and use_raw:
       x = x.raw.to_adata()
    
    if isinstance(x, MuData):
       x = x[mod]

    import pandas as pd
    out = pd.DataFrame({},  index=x.obs.index)
    
    coord = pd.DataFrame({})
    if obsm is not None:
        coord = pd.DataFrame(x.obsm[obsm], index=x.obs.index)
        coord = coord.rename(columns={0: "dim1", 1: "dim2"})
    
    meta = pd.DataFrame({})
    if obs is True:
      meta = pd.DataFrame(x.obs)
    elif isinstance(obs, str) or isinstance(obs, list):
      meta = pd.DataFrame(x.obs[obs])
    
    exprs = pd.DataFrame({})
    if features is not None:
      exprs = x[:, features].to_df(layer=layer)
    
    return out.join([coord, meta, exprs])