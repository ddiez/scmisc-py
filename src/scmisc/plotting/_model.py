def model_fit(model, train="train_loss_epoch", validation="validation_loss", ax=None, figsize=(4, 4), title=None):
  import matplotlib.pyplot as plt

  if ax is None:
    fig, ax = plt.subplots(figsize=figsize)
  
  ax.plot(model.history["train_loss_epoch"], label="train")
  ax.plot(model.history["validation_loss"], label="validation")
  if title is None:
    ax.set_title("Loss over training epochs")
  else:
    ax.set_title(title)
  ax.set_xlabel("Epoch")
  ax.set_ylabel("Loss")
  ax.legend()
  
  if ax is None:
    fig.show()