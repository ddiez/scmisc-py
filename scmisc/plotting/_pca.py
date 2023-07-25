import matplotlib.pyplot as plt
import numpy as np

def pca_variance_cumper(x):
  import numpy as np
  import matplotlib.pyplot as plt

  var = x.uns["pca"]["variance"]

  var = 100 * var / np.sum(var)
  var = np.cumsum(var)

  plt.scatter(x=range(var.shape[0]), y=var, c="black", s=1)
  plt.axhline(y=80, c="red", linewidth=1)
  plt.xlabel("Component")
  plt.ylabel("Cumulative percentage variance")

