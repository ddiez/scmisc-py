import matplotlib.pyplot as plt
from .._tools import calculate_dims

def volcano(x, groupby, N=10, lfc=1, nrows=None, ncols=None, figsize=None, lfc_col="lfc_mean", prob_col="proba_de", xlab="log2FC", ylab="-log10(p-value)"):
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
    tmp_up = tmp[tmp[lfc_col] > lfc]
    tmp_down = tmp[tmp[lfc_col] < -lfc]
    tmp_up = tmp_up.sort_values(prob_col, ascending=False).head(N)
    tmp_down = tmp_down.sort_values(prob_col, ascending=False).head(N)
    ax[crow, ccol].scatter(tmp[lfc_col], tmp[prob_col], s=1, c="black", marker=".")
    ax[crow, ccol].axvline(x=0, color="lightgray", linewidth=.5, linestyle="-")
    ax[crow, ccol].axvline(x=-1, color="black", linewidth=.5, linestyle="--")
    ax[crow, ccol].axvline(x=1, color="black", linewidth=.5, linestyle="--")
    ax[crow, ccol].set_title(group)
    ax[crow, ccol].set_axis_on()
    #ax[crow, ccol].set_xticks([])
    #ax[crow, ccol].set_yticks([])
    ax[crow, ccol].set_xlabel(xlab)
    ax[crow, ccol].set_ylabel(ylab)
    ax[crow, ccol].scatter(tmp_up[lfc_col], tmp_up[prob_col], s=1, c="red")
    ax[crow, ccol].scatter(tmp_down[lfc_col], tmp_down[prob_col], s=1, c="blue")
    text = []
    for i, d in tmp_up.iterrows():
      text.append(ax[crow, ccol].text(d[lfc_col], d[prob_col], s=i, c="red", size=8))

    for i, d in tmp_down.iterrows():
      text.append(ax[crow, ccol].text(d[lfc_col], d[prob_col], s=i, c="blue", size=8))

    adjust_text(text, arrowprops=dict(arrowstyle="-", color='k', lw=0.5), ax=ax[crow, ccol], precision=1, text_from_points=False)

    #ta.allocate_text(fig, ax[crow, ccol], x=tmp_up.lfc_mean, y=tmp_up.proba_de, text_list=tmp_up.index, textcolor="red", linecolor="black", min_distance=.2, max_distance=.3, textsize=8)
    #ta.allocate_text(fig, ax[crow, ccol], x=tmp_down.lfc_mean, y=tmp_down.proba_de, text_list=tmp_down.index, textcolor="blue", linecolor="black", min_distance=.2, max_distance=.3, textsize=8)

    ccol = ccol + 1
    if (ccol > ncols - 1):
      crow+=1
      ccol=0
  
  fig.tight_layout()

