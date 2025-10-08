def model_fit(model, train="train_loss_epoch", validation="validation_loss", ax=None, figsize=(4, 4), title=None):
  import matplotlib.pyplot as plt

  if ax is None:
    fig, ax = plt.subplots(figsize=figsize)

  history = model.history
  
  ax.plot(history[train], label=train)
  ax.plot(history[validation], label=validation)
  if title is None:
    ax.set_title("Loss over training epochs")
  else:
    ax.set_title(title)
  ax.set_xlabel("Epoch")
  ax.set_ylabel("Loss")
  ax.legend()
  
  if ax is None:
    fig.show()