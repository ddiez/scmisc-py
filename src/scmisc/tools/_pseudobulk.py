def pseudobulk(adata, groupby, splitby=None, aggregate="sum", layer=None, metadata=None, return_mudata=True):
  from anndata import AnnData
  from mudata import MuData
  import pandas as pd

  d = adata.to_df(layer=layer)
  d[groupby] = adata.obs[groupby]

  if splitby is not None:
    d[splitby] = adata.obs[splitby]
    d = d.groupby([groupby, splitby], observed=False)
  else:
    d = d.groupby(groupby, observed=False)

  match aggregate:
    case "sum":
      d = d.sum()
    case "mean":
      d = d.mean()
  
  groups = set(adata.obs[groupby].unique())

  obs = pd.DataFrame({"samplename": list(groups)}, index=list(groups))
  
  if splitby is not None:
    d = {k: v for k, v in d.groupby(splitby, observed=False)}

    for n in d:
      tmp = d[n].reset_index(splitby, drop=True)

      # Find missing groups and fill them with zeros.
      sel = list(set(groups).difference(tmp.index))
      for group in sel:
        tmp.loc[group, :] = 0

      d[n] = tmp
    
    if return_mudata:
      a = {}
      for n in d:
        a[n] = AnnData(d[n].set_index(obs.index), obs=obs, var=adata.var)
      
      res = MuData(a)
    else: res = d

  else:
    res = MuData({
      "pseudobulk": AnnData(d, obs=obs, var=adata.var)
    })

  if metadata is not None:
    for n in res.mod_names:
      res[n].obs = metadata.loc[res[n].obs_names]
  
  return res