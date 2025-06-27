import matplotlib.pyplot as plt
import numpy as np

def pca_variance_cumper(x):
  import numpy as np
  import matplotlib.pyplot as plt

  var = x.uns["pca"]["variance"]

  var = 100 * var / np.sum(var)
  var = np.cumsum(var)

  fig, ax = plt.subplots()
  ax.scatter(x=range(var.shape[0]), y=var, c="black", s=1)
  ax.axhline(y=80, c="red", linewidth=1)
  ax.set_xlabel("Component")
  ax.set_ylabel("Cumulative percentage variance")
  ax.xaxis.grid(True, color="gray", linestyle="--", linewidth=0.5)

