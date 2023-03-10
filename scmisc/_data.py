import pkg_resources

def load_10x_doublet_rate():
  import pandas as pd
  stream = pkg_resources.resource_stream(__name__, 'datasets/10x_doublets.csv')
  return pd.read_csv(stream)