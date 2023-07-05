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

def plot_xy(data, x, y, expand=None, size=.1, color="lightgrey", highlight_color="red", nrows=None, ncols=None, figsize=None, grid=False, axis=False, return_ax=False, *args, **kwargs):
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

    fig.tight_layout()
  
  if return_ax:
    return ax

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

def plot_volcano(x, groupby, N=10, lfc=1, nrows=None, ncols=None, figsize=None):
  from adjustText import adjust_text
  #import textalloc as ta

  groups = x[groupby].unique()
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
    tmp = x[x[groupby] == group]
    tmp_up = tmp[tmp.lfc_mean > lfc]
    tmp_down = tmp[tmp.lfc_mean < -lfc]
    tmp_up = tmp_up.sort_values("proba_de", ascending=False).head(N)
    tmp_down = tmp_down.sort_values("proba_de", ascending=False).head(N)
    ax[crow, ccol].scatter(tmp.lfc_mean, tmp.proba_de, s=1, c="black", marker=".")
    ax[crow, ccol].axvline(x=0, color="lightgray", linewidth=.5, linestyle="-")
    ax[crow, ccol].axvline(x=-1, color="black", linewidth=.5, linestyle="--")
    ax[crow, ccol].axvline(x=1, color="black", linewidth=.5, linestyle="--")
    ax[crow, ccol].set_title(group)
    ax[crow, ccol].set_axis_on()
    ax[crow, ccol].set_xticks([])
    ax[crow, ccol].set_yticks([])
    ax[crow, ccol].set_xlabel("log10 FC")
    ax[crow, ccol].set_ylabel("Probability DE")
    ax[crow, ccol].scatter(tmp_up.lfc_mean, tmp_up.proba_de, s=1, c="red")
    ax[crow, ccol].scatter(tmp_down.lfc_mean, tmp_down.proba_de, s=1, c="blue")
    text = []
    for i, d in tmp_up.iterrows():
      text.append(ax[crow, ccol].text(d.lfc_mean, d.proba_de, s=i, c="red", size=8))

    for i, d in tmp_down.iterrows():
      text.append(ax[crow, ccol].text(d.lfc_mean, d.proba_de, s=i, c="blue", size=8))

    adjust_text(text, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), ax=ax[crow, ccol], precision=1, text_from_points=False)

    #ta.allocate_text(fig, ax[crow, ccol], x=tmp_up.lfc_mean, y=tmp_up.proba_de, text_list=tmp_up.index, textcolor="red", linecolor="black", min_distance=.2, max_distance=.3, textsize=8)
    #ta.allocate_text(fig, ax[crow, ccol], x=tmp_down.lfc_mean, y=tmp_down.proba_de, text_list=tmp_down.index, textcolor="blue", linecolor="black", min_distance=.2, max_distance=.3, textsize=8)

    ccol = ccol + 1
    if (ccol > ncols - 1):
      crow+=1
      ccol=0
  
  fig.tight_layout()