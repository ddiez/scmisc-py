import matplotlib.pyplot as plt
import numpy as np

def calculate_dims(ngroups, nrows=None, ncols=None):
  if ncols is None and nrows is None:
      ncols = np.int32(np.ceil(np.sqrt(ngroups)))
      nrows = int(np.ceil(ngroups/ncols))

  if ncols is None and nrows is not None:
    ncols = int(np.ceil(ngroups/nrows))

  if nrows is None and ncols is not None:
    nrows = int(np.ceil(ngroups/ncols))

  return nrows, ncols

def xy(data, x, y, expand=None, size=.1, color="lightgrey", highlight_color="red", nrows=None, ncols=None, figsize=None, grid=False, axis=False, return_ax=False, *args, **kwargs):
  X = data[x]
  Y = data[y]

  if expand is None:
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, squeeze=True)
    ax.scatter(X, Y, s=size, c=color, *args, **kwargs)
  else:
    expand = data[expand]
    groups = expand.unique()
    ngroups = groups.shape[0]

    nrows, ncols = calculate_dims(ngroups=ngroups, nrows=nrows, ncols=ncols)
    
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
      cells = index
      
      ax[crow, ccol].scatter(X, Y, s=size, c=color, marker="o")
      ax[crow, ccol].scatter(X[cells], Y[cells], s=size, c=highlight_color, marker="o")
      ax[crow, ccol].set_title(str(groups[k]))

      if axis:
        ax[crow, ccol].set_axis_on()

      if grid:
        ax[crow, ccol].grid(True)
      
      ccol = ccol + 1
      if (ccol > ncols - 1):
        crow+=1
        ccol=0

    fig.supxlabel(x)
    fig.supylabel(y)
    fig.tight_layout()
  
  if return_ax:
    return ax

def coord(x, expand=None, basis=None, size=.1, color="lightgrey", highlight_color="red", nrows=None, ncols=None, figsize=None):  
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

    nrows, ncols = calculate_dims(ngroups=ngroups, nrows=nrows, ncols=ncols)

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
