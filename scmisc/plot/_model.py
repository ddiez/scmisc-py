import matplotlib.pyplot as plt
import numpy as np

def plot_model_fit(model, train="train_loss_epoch", validation="validation_loss"):
  plt.plot(model.history["train_loss_epoch"], label="train")
  plt.plot(model.history["validation_loss"], label="validation")
  plt.title("Loss over training epochs")
  plt.xlabel("Epoch")
  plt.ylabel("Loss")
  plt.legend()
  plt.show()