def pseudobulk(adata, groupby, splitby=None, aggregate="mean"):
  from anndata import AnnData
  from mudata import MuData

  d = adata.to_df()
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
  
  if splitby is not None:
    d = {k: v for k, v in d.groupby(splitby, observed=False)}

    for n in d:
      tmp = d[n].reset_index(splitby)
      tmp = tmp.drop(columns=[splitby])
      d[n] = tmp
    
    a = {}
    for n in d:
      a[n] = AnnData(d[n])
    res = MuData(a)
    
  else:
    res = AnnData(d)
  
  return res