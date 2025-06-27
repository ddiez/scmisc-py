def fetch_data(x, obsm=None, obs=None, features=None, use_raw=None):
    if use_raw is True:
       x = x.raw.to_adata()

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
      exprs = x[:, features].to_df()
    
    return out.join([coord, meta, exprs])