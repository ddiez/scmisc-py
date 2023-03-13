def predict_doublet_rate(x, batch=None, fit=None):
  import numpy as np
  import pandas as pd

  if fit is None:
      fit = generate_doublet_rate_model_10x()

  import anndata
  if isinstance(x, anndata.AnnData):
    x = x.obs

  res = []

  if batch is None:
    N = x.shape[0]
    n = np.array(N).reshape(-1, 1)
    res.append({"batch": batch, "count": N, "doublet_rate": fit.predict(n)[0, 0]})

  if batch is not None:
    for b in x[batch].values.categories:
      N = x[x[batch] == b].shape[0]
      n = np.array(N).reshape(-1, 1)
      rate = fit.predict(n)[0, 0]
      n_doublet = N * rate / 100
      tmp = {"batch": b, "count": N, "doublet_rate": rate, "n_doublet": n_doublet}
      res.append(tmp)
  
  return pd.DataFrame(res)

def generate_doublet_rate_model_10x(data=None):
  if data is None:
      from ._data import load_10x_doublet_rate
      data = load_10x_doublet_rate()
  
  from sklearn.linear_model import LinearRegression
  
  X = data["count"].to_numpy().reshape(-1,1)
  y = data["doublet_rate"].to_numpy().reshape(-1,1)
  fit = LinearRegression().fit(X, y)
  return fit

def plot_doublet_ratio(x=None):
  from ._data import load_10x_doublet_rate
  data = load_10x_doublet_rate()
  fit = generate_doublet_rate_model_10x(data)

  ax = data.plot("count", "doublet_rate", kind="scatter", color="black")
  ax.axline((0, fit.intercept_[0]), slope=fit.coef_[0,0], color="red", linewidth=1)
  if x is not None:
    ax.scatter(x["count"], x["doublet_rate"], marker="x", color="blue")
