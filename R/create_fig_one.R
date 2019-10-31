#' Creates Simulation Figure 1
#'
#' This function creates a table comparing the novel NDP
#' to a classic heirarchical finite mixture model
#'
#' @param K_f *vector* of the different number of cluster components in fhmm
#' @param L_f number of within-cluster mixture components in fhmm
#' @param K   number of cluster sticks for NDP
#' @param L number of within-cluster mixture-sticks for NDP
#' @param J total number of groups
#' @param f_s list of functions that should be the same length as K_true
#' @param iter_max maximum number of iterations to use in model fitting
#' @param warm_up number of warm-up iterations to use in model fitting
#' @return plot showing results
#' @export
create_fig_one <- function(num_sims = 5,
                             K_f = 2:8,
                             L_f = 4,
                             K = 10,
                             L = 4,
                             J = c(50,50,50),
                             R = 5,
                             f_s = c(function(x) 3*((1/2)*dbeta(x/R,1,8) + (1/2)*dbeta(x/R,6,1)),
                                     function(x) 3*((1/5)*dbeta(x/R,3,2) + (2/3)*dbeta(x/R,3,1) + (2/15)*dbeta(x/R,1,1)),
                                     function(x) 3*((1/2)*dbeta(x/R,8,2) + (1/2)*dbeta(x/R,30,50))),
                             iter_max = 1.5E4,
                             warm_up = 1.2E4,
                             thin = 1,
                             seed = NULL) {

  # Preliminary Checks ------------------------------------------------------

  stopifnot(length(f_s) == length(J))
  if(is.null(seed))
    set.seed(341351)
  else
    set.seed(seed)

  pltdf <- purrr::map_dfr(K_f,function(x) create_table_one(K_f = x,K=8)$raw_table %>%
                            dplyr::mutate(K = x))

  p <- pltdf %>%
    tidyr::gather(loss,sqr_err,key="Metric",value="Metric_Value") %>%
    dplyr::group_by(K,model,Metric) %>%
    dplyr::summarise(Mean_value = mean(Metric_Value),
              sd_value = sd(Metric_Value)/sqrt(dplyr::n())) %>%
    dplyr::ungroup() %>% dplyr::mutate(upper = Mean_value + 2 * sd_value,
                                       lower = Mean_value - 2 * sd_value) %>%
    ggplot2::ggplot(ggplot2::aes(x=K,y=Mean_value,fill=model)) +
    ggplot2::geom_point() + ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=0.35) + ggplot2::facet_wrap(~Metric,scales = "free_y") +
    ggplot2::theme_bw() + ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(breaks=K_f) +
    ggplot2::labs(title = "Comparison of Loss Metrics Between Models, Across Number of Cluster Functions",
                  y = "Error")


  return(p)

}
