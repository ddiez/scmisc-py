def plot_coord(x, expand=None, basis=None, size=.1, color="lightgrey", highlight_color="red", ncols=None, figsize=None):
  import matplotlib.pyplot as plt
  import numpy as np
  
  if basis is None:
    basis = adata.obsm_keys()[0]
  else:
    if basis not in x.obsm_keys():
      basis = f"X_{basis}"
    if basis not in x.obsm_keys():
      return(f"Error, basis {basis} not found.")
    
  key = basis.removeprefix("X_").upper()
  keys = (f"{key}1", f"{key}2")
  
  coord = x.obsm[basis]
  
  if expand is not None:
    expand = x.obs[expand].sort_values()
    groups = expand.unique()
    ngroups = groups.shape[0]
    
    if ncols is None:
      ncols = np.int32(np.ceil(np.sqrt(ngroups)))
    
    nrows = int(np.ceil(ngroups/ncols))
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in range(nrows):
      for j in range(ncols):
        ax[i, j].grid(False)
        ax[i, j].set_axis_off()

    crow=0
    ccol=0
    
    for k in range(ngroups):
      group = groups[k]
      index = expand[expand == group].index
      cells = x.obs_names.isin(index)
      
      im = ax[crow, ccol].scatter(coord[:,0], coord[:,1], s=size, c=color, marker="o")
      ax[crow, ccol].scatter(coord[cells, 0], coord[cells, 1], s=size, c=highlight_color, marker="o")
      ax[crow, ccol].set_title(str(groups[k]))
      ax[crow, ccol].set_axis_on()
      ax[crow, ccol].set_xticks([])
      ax[crow, ccol].set_yticks([])
      ax[crow, ccol].set_xlabel(keys[0])
      ax[crow, ccol].set_ylabel(keys[1])
      
      ccol = ccol + 1
      if (ccol > ncols - 1):
        crow+=1
        ccol=0

    fig.tight_layout()
  else:
    plt.scatter(coord[:, 0], coord[:, 1], s=size, c=color, marker="o")

  plt.show()