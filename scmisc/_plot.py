import matplotlib.pyplot as plt
import numpy as np

def plot_xy(x, y, expand=None, size=.1, color="lightgrey", highlight_color="red", nrows=None, ncols=None, figsize=None, *args, **kwargs):
  if expand is None:
    plt.scatter(x, y, *args, **kwargs)
  else:
    groups = expand.unique()
    ngroups = groups.shape[0]

    if ncols is None and nrows is None:
      ncols = np.int32(np.ceil(np.sqrt(ngroups)))
      nrows = int(np.ceil(ngroups/ncols))

    if ncols is None and nrows is not None:
      ncols = int(np.ceil(ngroups/nrows))

    if nrows is None and ncols is not None:
      nrows = int(np.ceil(ngroups/ncols))
    
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, squeeze=False)

    for i in range(nrows):
      for j in range(ncols):
        ax[i, j].grid(False)
        ax[i, j].set_axis_off()

    crow=0
    ccol=0
    
    for k in range(ngroups):
      group = groups[k]
      index = expand[expand == group].index
      #cells = x.obs_names.isin(index)
      cells = index
      
      ax[crow, ccol].scatter(x, y, s=size, c=color, marker="o")
      ax[crow, ccol].scatter(x[cells], y[cells], s=size, c=highlight_color, marker="o")
      ax[crow, ccol].set_title(str(groups[k]))
      ax[crow, ccol].set_axis_on()
      ax[crow, ccol].set_xticks([])
      ax[crow, ccol].set_yticks([])
      #ax[crow, ccol].set_xlabel(keys[0])
      #ax[crow, ccol].set_ylabel(keys[1])
      
      ccol = ccol + 1
      if (ccol > ncols - 1):
        crow+=1
        ccol=0

    fig.tight_layout()

def plot_coord(x, expand=None, basis=None, size=.1, color="lightgrey", highlight_color="red", nrows=None, ncols=None, figsize=None):  
  if basis is None:
    basis = x.obsm_keys()[0]
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
    
    if ncols is None and nrows is None:
      ncols = np.int32(np.ceil(np.sqrt(ngroups)))
      nrows = int(np.ceil(ngroups/ncols))

    if ncols is None and nrows is not None:
      ncols = int(np.ceil(ngroups/nrows))

    if nrows is None and ncols is not None:
      nrows = int(np.ceil(ngroups/ncols))
    
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, squeeze=False)
    
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
      
      ax[crow, ccol].scatter(coord[:,0], coord[:,1], s=size, c=color, marker="o")
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

def plot_pca_variance_cumper(x):
  import numpy as np
  import matplotlib.pyplot as plt

  var = x.uns["pca"]["variance"]

  var = 100 * var / np.sum(var)
  var = np.cumsum(var)

  plt.scatter(x=range(var.shape[0]), y=var, c="black", s=1)
  plt.axhline(y=80, c="red", linewidth=1)
  plt.xlabel("Component")
  plt.ylabel("Cumulative percentage variance")
