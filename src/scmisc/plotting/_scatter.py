import matplotlib.pyplot as plt
from pandas import CategoricalDtype
from .._tools import calculate_dims
from numpy import ndarray

def xy(data, x, y, expand=None, size=1, color_density=False, color="lightgrey", highlight_color="red", nrows=None, ncols=None, figsize=None, grid=False, axis=False, return_ax=False, *args, **kwargs):
  X = data[x]
  Y = data[y]

  if color_density:
    from sklearn.neighbors import KernelDensity
    kd = KernelDensity()

    fit = kd.fit(data[[x, y]])
    data[".density"] = fit.score_samples(data[[x, y]])
    color = data[".density"]

  if expand is None:
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, squeeze=True)
    ax.scatter(X, Y, s=size, c=color, *args, **kwargs)
  else:
    expand = data[expand]

    if isinstance(expand.dtype, CategoricalDtype):
      groups = expand.cat.categories
    else:
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
      
      ax[crow, ccol].scatter(X, Y, s=size, c=color, marker=".", linewidth=0)
      ax[crow, ccol].scatter(X[cells], Y[cells], s=size, c=highlight_color, marker=".", linewidth=0)
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

def coord(x, expand=None, groups=None, basis=None, size=1, color="lightgrey", highlight_color="red", nrows=None, ncols=None, figsize=None):  
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

    if groups is None:
      groups = expand.unique()

    if isinstance(groups, list):
      ngroups = len(groups)
    
    if isinstance(groups, ndarray):
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
      
      ax[crow, ccol].scatter(coord[:,0], coord[:,1], s=size, c=color, marker=".", linewidth=0)
      ax[crow, ccol].scatter(coord[cells, 0], coord[cells, 1], s=size, c=highlight_color, marker=".", linewidth=0)
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
    plt.scatter(coord[:, 0], coord[:, 1], s=size, c=color, marker=".", linewidth=0)

  plt.show()
