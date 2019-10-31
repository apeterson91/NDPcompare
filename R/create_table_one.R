#' Simulation Table 1
#'
#' This function creates a table comparing the novel NDP
#' to a classic heirarchical finite mixture model
#'
#' @param K_f number of cluster components in fhmm
#' @param L_f number of within-cluster mixture components in fhmm
#' @param K number of cluster sticks for NDP
#' @param L number of within-cluster mixture-sticks for NDP
#' @param J total number of groups
#' @param f_s list of functions that should be the same length as K_true
#' @param iter_max maximum number of iterations to use in model fitting
#' @param warm_up number of warm-up iterations to use in model fitting
#' @return list containing raw and aggregated table values
#' @export
create_table_one <- function(num_sims = 5,
                             K_f = 4,
                             L_f = 4,
                             K = 4,
                             L = 4,
                             J = c(50,50,50),
                             R = 5,
                             f_s = c(function(x) 3*((1/2)*dbeta(x/R,1,8) + (1/2)*dbeta(x/R,6,1)),
                                      function(x) 3*((1/5)*dbeta(x/R,3,2) + (2/3)*dbeta(x/R,3,1) + (2/15)*dbeta(x/R,1,1)),
                                      function(x) 3*((1/2)*dbeta(x/R,8,2) + (1/2)*dbeta(x/R,30,50))),
                             iter_max = 1.5E4,
                             warm_up = 1.3E4,
                             thin = 1,
                             seed = NULL) {

  # Preliminary Checks ------------------------------------------------------

  stopifnot(length(f_s) == length(J))
  if(is.null(seed))
    set.seed(341351)
  else
    set.seed(seed)

  # Generate True Adjacency Matrix ------------------------------------------

  A <- generate_adjacency_matrix(J)

  # Simulate Data -----------------------------------------------------------

  data <- purrr::map(1:num_sims,function(y){

    df <- purrr::map_dfr(1:length(J), function(x) rndpp::rnhpp(nsim = J[x],
                                                               lambda = f_s[[x]],
                                                               interval = c(0,R),
                                                               max=max(f_s[[x]](seq(from=0,to=R,by=0.01)) ) ) %>%
                           dplyr::mutate(distance_id = sim_id + J[x]*(x-1),
                                         group_id = x,
                                         simulation_id = y))
  })

  ndp_fits <- purrr::map2_dfr(1:length(data),data,function(x,y){
    r <- y %>% dplyr::arrange(distance_id) %>%
      dplyr::select(event_times) %>% dplyr::pull()

    n_j <- y %>% dplyr::arrange(distance_id) %>%
      dplyr::group_by(distance_id) %>% dplyr::count() %>%
      dplyr::ungroup() %>% dplyr::mutate(start = (cumsum(n) ) ) %>%
      dplyr::mutate(start_ = tidyr::replace_na(dplyr::lag(start),0) ) %>% dplyr::select(-start) %>%
      dplyr::rename(start=start_,go =n) %>%
      dplyr::select(start,go) %>% as.matrix()

    NDP_fit <- rndpp::nd_nhpp_fixed(X = as.matrix(rep(1,nrow(n_j))),
                                    n_j = n_j,r = r,
                                    mu_0 = 0,kappa_0 = 1,nu_0 = 1,
                                    sigma_0 = 1,
                                    L = L,
                                    K = K,
                                    alpha = 1,
                                    rho = 1,
                                    iter_max = iter_max,
                                    warm_up = warm_up,
                                    thin = thin,
                                    seed = 34134+x)

    out <- dplyr::tibble(sim = x,
                         loss = rndpp::green_loss(NDP_fit,truth=A)$loss,
                         sqr_err = sum((A-NDP_fit$pmat)^2),
                         model = "NDP")

  })

  fhm_fits <- purrr::map2_dfr(1:length(data),data,function(x,y){
    r <- y %>% dplyr::arrange(distance_id) %>%
      dplyr::select(event_times) %>% dplyr::pull()

    n_j <- y %>% dplyr::arrange(distance_id) %>%
      dplyr::group_by(distance_id) %>% dplyr::count() %>%
      dplyr::ungroup() %>% dplyr::mutate(start = (cumsum(n) ) ) %>%
      dplyr::mutate(start_ = tidyr::replace_na(dplyr::lag(start),0) ) %>% dplyr::select(-start) %>%
      dplyr::rename(start=start_,go =n) %>%
      dplyr::select(start,go) %>% as.matrix()

    fhm_fit <- fhmm::fhmm(n_j = n_j,r = r,
                          mu_0 = 0,
                          kappa_0 = 1,nu_0 = 1,
                          sigma_0 = 1,
                          L = L_f,
                          K = K_f,
                          iter_max = iter_max,
                          warm_up = warm_up,
                          thin = thin,
                          seed = 34134)

    out <- dplyr::tibble(sim = x,
                         loss = fhmm::green_loss(fhm_fit,truth = A)$loss,
                         sqr_err = sum((A-fhm_fit$pmat)^2),
                         model = "FHMM")

  })


  raw_table <- rbind(ndp_fits,fhm_fits)
  raw_table <- raw_table %>% dplyr::mutate(model = factor(model))

  require(tables)
  agg_table <- tabular( (loss + sqr_err)*(mean+sd) ~ (model),data=raw_table)

  out <- list(raw_table = raw_table,
              agg_table = agg_table)

  return(out)

}



