def filter_cells_by_genes(adata, genes, cutoff=0, layer="counts"):
  for gene in genes:
    if gene in adata.var_names:
      sel = adata[:, adata.var_names == gene].layers[layer].toarray() > cutoff
      adata = adata[~sel, :]
  return adata

def calculate_dims(ngroups, nrows=None, ncols=None):
  import numpy as np

  if ncols is None and nrows is None:
      ncols = np.int32(np.ceil(np.sqrt(ngroups)))
      nrows = np.int32(np.ceil(ngroups/ncols))

  if ncols is None and nrows is not None:
    ncols = np.int32(np.ceil(ngroups/nrows))

  if nrows is None and ncols is not None:
    nrows = np.int32(np.ceil(ngroups/ncols))

  return nrows, ncols

def get_next_ax(row, col, ncols):
  ''' 
    # Get the next row and column for a grid of subplots
  '''
  col = col + 1
  if (col > ncols - 1):
    row+=1
    col=0
  
  return row, col