import scipy.sparse

class sparse_crs(scipy.sparse.csr_array):
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    
  def sum(self, axis=None, keepdims=False, **kwargs):
    s = super().sum(axis=axis, **kwargs)
    
    if keepdims:
      return s.reshape(-1, 1)
    else:
      return s