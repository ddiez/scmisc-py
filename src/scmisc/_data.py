
def load_10x_doublet_rate():
  from importlib.resources import files
  import pandas as pd
  
  path = files(__name__).joinpath("datasets/10x_doublets.csv")
  return pd.read_csv(path)