###
###
###
###    Purpose:   Top-level functions to compute economic values
###    started:   2018-10-29 (pvr)
###
### ############################################################## ###


#' @title Compute Economic Values Based On Revenue
#'
#' @description
#' In this function, the economic value is computed as the change
#' in revenue for a given change in the population mean. The computation
#' of the economic value is simplified by the assumption that the
#' costs can be kept constant for the different levels of the population
#' mean.
#'
#' @param pn_mean         empirical population mean
#' @param pn_sd           empirical population standard deviation
#' @param pvec_class_freq vector of empirical trait frequencies
#' @param pvec_threshold  vector of thresholds given by payment system
#' @param pvec_price      vector of prices given by payment system
#' @param pn_delta_mean   small change of population mean
#'
#' @return
#' @export
compute_economic_value <- function(pn_mean,
                                   pn_sd,
                                   pvec_class_freq = NULL,
                                   pvec_threshold  = NULL,
                                   pvec_price,
                                   pn_delta_mean  ){
  ### # checking input parameters
  if ( length(pn_mean) != 1L )
    stop("Population mean (pn_mean) must be a single number: ", length(pn_mean))
  if ( length(pn_sd) != 1L )
    stop("Population standard deviation (pn_sd) must be a single number: ", length(pn_sd))
  if ( length(pn_delta_mean) != 1L )
    stop("Change in population mean (pn_delta_mean) must be a single number: ", length(pn_delta_mean))

  ### # for a continuous trait pvec_class_freq is NULL
  if ( is.null(pvec_class_freq) ) {
    ### # for continuous distributions, we need a vector of thresholds
    if ( is.null(pvec_threshold) ) {
      stop("For a continuous trait, vector of thresholds for price categories cannot be NULL")
    }
    ### # vector of prices must contain one element more than vector of thresholds
    if ( length(pvec_threshold) != (length(pvec_price)-1) )
      stop("Number of thresholds and number prices not consistent: ", length(pvec_threshold), " ", length(pvec_price))

    ### # use the specified thresholds for computing economic values later
    vec_threshold <- pvec_threshold

  } else {
    ### # case of a discrete trait, check that number of class frequencies and number of prices are the same
    if ( length(pvec_class_freq) != length(pvec_price) )
      stop("Number of classes must be the same as number of prices: ", length(pvec_class_freq), " ", length(pvec_price))

    ### # thresholds assuming normal distribution, if any thresholds are specified, they are ignored
    vec_cum_freq_dist <- cumsum(pvec_class_freq)
    vec_threshold <- qnorm(vec_cum_freq_dist[1:(length(vec_cum_freq_dist)-1)], mean = pn_mean, sd = pn_sd, lower.tail = TRUE)

  }

  ### # compute economic value based on mean, sd, threshold, prices and delta
  ev_result <- compute_ev_base_rev( pn_mean        = pn_mean,
                                    pn_sd          = pn_sd,
                                    pvec_threshold = vec_threshold,
                                    pvec_price     = pvec_price,
                                    pn_delta_mean  = pn_delta_mean  )

  return(ev_result)

}


#' @title Compute numeric value of economic value based on revenue
#'
#' @description
#' This function computes the economic value based on changes in revenue for
#' a small change in the population mean.
#'
#' @param pn_mean         empirical population mean
#' @param pn_sd           empirical population standard deviation
#' @param pvec_threshold  vector of thresholds given by payment system
#' @param pvec_price      vector of prices given by payment system
#' @param pn_delta_mean   small change of population mean
#' @return ev_rev_result
compute_ev_base_rev <- function( pn_mean,
                                 pn_sd,
                                 pvec_threshold,
                                 pvec_price,
                                 pn_delta_mean ) {
  ### # getting area under normal distribution for base situation
  vec_freq_base <- pnorm(pvec_threshold, mean = pn_mean, sd = pn_sd)
  vec_freq_base <- c(vec_freq_base, 1-sum(vec_freq_base))
  ### # compute the expected price in the base situation
  n_ex_price_base <- crossprod(vec_freq_base, pvec_price)

  ### # shift the distribution
  n_mean_shifted <- pn_mean + pn_delta_mean
  vec_freq_shifted <- pnorm(pvec_threshold, mean = n_mean_shifted, sd = pn_sd)
  vec_freq_shifted <- c(vec_freq_shifted, 1-sum(vec_freq_shifted))
  n_ex_price_shifted <- crossprod(vec_freq_shifted, vec_price_cca)

  ### # difference corresponds to ev
  ev_result <- n_ex_price_shifted - n_ex_price_base

  return(ev_result)
}
