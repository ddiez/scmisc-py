def purity(adata, x, y):
  from plotnine import ggplot, aes, geom_tile, scale_fill_gradient, theme, element_text

  obs = adata.obs[[x, y]]
  obs = obs.groupby(y).value_counts(normalize=True).reset_index(name="purity")

  return (ggplot(obs, aes(x, y, fill="purity")) + 
   geom_tile() +
   scale_fill_gradient(low="white", high="red") +
   theme(figure_size=(2, 5),
         axis_text=element_text(color="black"),
         axis_text_x=element_text(rotation=45, hjust=1, vjust=1)))
