##################
# cor_within_env #
##################

cor_within_env <- function(d, obs = "trait", pred = "pred") {
  # within environment correlation
  cor_with_env <- d %>%
    group_by(Env) %>%
    summarise(cor = cor(.data[[obs]], .data[[pred]], use = "pairwise.complete.obs"))
  
  return(cor_with_env)
}