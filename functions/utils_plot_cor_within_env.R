#######################
# plot_cor_within_env #
#######################

plot_cor_within_env <- function(d, obs = "trait", pred = "pred", title) {
  
  p <- ggplot(data = d, aes(x = pred, y = trait)) +
    geom_point() +
    facet_wrap(~Env, nrow = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    ggtitle(title)
  
  return(p)
}