def rank_vector(x, descending=True):
  tmp = np.argsort(x)
  if descending:
    tmp = np.flip(tmp)
  return(np.argsort(tmp))
