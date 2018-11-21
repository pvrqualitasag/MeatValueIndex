###
###
###
###    Purpose:   Top-level functions to compute economic values
###    started:   2018-10-29 (pvr)
###
### ############################################################## ###


#' @title Compute Economic Values Based On Price Changes
#'
#' @description
#' In this function, the economic value is computed as the change
#' in expected average price over a population for a given change
#' in the population mean. The computation
#' of the economic value is simplified by the assumption that the
#' costs can be kept constant for the different levels of the population
#' mean. Futhermore, we assume that the change in population mean
#' does not change the amount of product produced. The result is
#' standardized and expressed per unit of the trait.
#'
#' @param pn_mean         empirical population mean
#' @param pn_sd           empirical population standard deviation
#' @param pvec_class_freq vector of empirical trait frequencies
#' @param pvec_threshold  vector of thresholds given by payment system
#' @param pvec_price      vector of prices given by payment system
#' @param pn_delta_mean   small change of population mean
#' @param pn_gen_sd       genetic standard deviation of trait
#' @param pb_verbose      level of output verbosity
#'
#' @return ev_result      list of two economic values, per unit trait and per genetic sd, if pn_gen_sd is given
#' @export compute_economic_value
compute_economic_value <- function(pn_mean,
                                   pn_sd,
                                   pvec_class_freq = NULL,
                                   pvec_threshold  = NULL,
                                   pvec_price,
                                   pn_delta_mean,
                                   pn_gen_sd       = NULL,
                                   pb_verbose      = FALSE ){
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

    if (pb_verbose) {
      cat("[INFO -- compute_economic_value] Continuous trait: threshold:\n")
      print(vec_threshold)
    }

    ### # compute economic value based on mean, sd, threshold, prices and delta
    ev_result <- compute_ev_price_cont( pn_mean        = pn_mean,
                                        pn_sd          = pn_sd,
                                        pvec_threshold = vec_threshold,
                                        pvec_price     = pvec_price,
                                        pn_delta_mean  = pn_delta_mean,
                                        pn_gen_sd      = pn_gen_sd,
                                        pb_verbose     = pb_verbose )

  } else {
    ### # case of a discrete trait, check that number of class frequencies and number of prices are the same
    if ( length(pvec_class_freq) != length(pvec_price) )
      stop("Number of classes must be the same as number of prices: ", length(pvec_class_freq), " ", length(pvec_price))

    ### # thresholds assuming normal distribution, if any thresholds are specified, they are ignored
    vec_cum_freq_dist <- cumsum(pvec_class_freq)
    vec_threshold <- qnorm(vec_cum_freq_dist[1:(length(vec_cum_freq_dist)-1)], mean = pn_mean, sd = pn_sd, lower.tail = TRUE)

    if (pb_verbose) {
      cat("[INFO -- compute_economic_value] Class freq:\n")
      print(pvec_class_freq)
      cat("[INFO -- compute_economic_value] price:\n")
      print(pvec_price)
      cat("[INFO -- compute_economic_value] Discrete trait: threshold:\n")
      print(vec_threshold)
      cat("[INFO -- compute_economic_value] Cumulative frequencies: \n")
      print(vec_cum_freq_dist)
    }

    ### # compute economic value based on discrete distribution given by pvec_class_freq
    ev_result <- compute_ev_price_disc( pn_mean         = pn_mean,
                                        pn_sd           = pn_sd,
                                        pvec_class_freq = pvec_class_freq,
                                        pvec_threshold  = vec_threshold,
                                        pvec_price      = pvec_price,
                                        pn_delta_mean   = pn_delta_mean,
                                        pn_gen_sd       = pn_gen_sd,
                                        pb_verbose      = pb_verbose  )
  }

  ### # return result
  return(ev_result)

}


#' @title Compute economic value based on price change for discrete trait
#'
#' @description
#' Economic value is computed based on price change between base distribution
#' and a shifted distribution for a discrete trait.
#'
#' @param pn_mean         empirical population mean
#' @param pn_sd           empirical population standard deviation
#' @param pvec_class_freq vector of empirical trait frequencies
#' @param pvec_threshold  vector of thresholds given by payment system
#' @param pvec_price      vector of prices given by payment system
#' @param pn_delta_mean   small change of population mean
#' @param pn_gen_sd       genetic standard deviation of trait
#' @param pb_verbose      level of output verbosity
#' @return ev_result      list of two economic values, per unit trait and per genetic sd, if pn_gen_sd is given
compute_ev_price_disc <- function( pn_mean,
                                   pn_sd,
                                   pvec_class_freq,
                                   pvec_threshold,
                                   pvec_price,
                                   pn_delta_mean,
                                   pn_gen_sd,
                                   pb_verbose       = FALSE ){
  ### # compute the expected price in the base situation
  n_ex_price_base <- crossprod(pvec_class_freq, pvec_price)
  if (pb_verbose) cat("[INFO -- compute_ev_price_disc] Base price: ", n_ex_price_base, "\n")
  ### # shift the distribution
  n_mean_shifted <- pn_mean + pn_delta_mean
  if (pb_verbose) cat("[INFO -- compute_ev_price_disc] Shifted mean: ", n_mean_shifted, "\n")
  ### # get cumulative frequencies under shifted distribution
  vec_freq_shifted_cum <- pnorm(pvec_threshold, mean = n_mean_shifted, sd = pn_sd)
  vec_freq_shifted <- diff(c(0, vec_freq_shifted_cum))
  if (pb_verbose) {
    cat("[INFO -- compute_ev_price_disc] Shifted freq cum\n")
    print(vec_freq_shifted_cum)
    cat("[INFO -- compute_ev_price_disc] Shifted freq diff\n")
    print(vec_freq_shifted)
  }
  ### # extract the sum of the shifted frequencies as last element of cumulative frequencies
  n_sum_freq_shifted <- vec_freq_shifted_cum[length(vec_freq_shifted_cum)]
  ### # add last entry for frequency
  if (n_sum_freq_shifted > 1){
    stop(" *** Error: compute_ev_price_disc: sum of shifted frequencies > 1: ", n_sum_freq_shifted)
  } else {
    vec_freq_shifted <- c(vec_freq_shifted, 1-n_sum_freq_shifted)
  }
  if (pb_verbose) {
    cat("[INFO -- compute_ev_price_disc] Shifted freq diff\n")
    print(vec_freq_shifted)
  }
  ### # compute price under shifted distribution
  n_ex_price_shifted <- crossprod(vec_freq_shifted, pvec_price)
  if (pb_verbose) cat("[INFO -- compute_ev_price_disc] Shifted price: ", n_ex_price_shifted, "\n")

  ### # difference corresponds to ev per trait unit, when rescaled with pn_delta_mean
  ev_result_per_trait_unit <- (n_ex_price_shifted - n_ex_price_base) / pn_delta_mean
  ### # if genetic standard deviation is specified, compute the ev on the basis of one genetic sd
  ev_result_per_gen_sd <- NULL
  if (!is.null(pn_gen_sd))
    ev_result_per_gen_sd <- ev_result_per_trait_unit * pn_gen_sd
  ### # return the results as a list
  ev_result <- list(ev_per_trait_unit = ev_result_per_trait_unit,
                    ev_per_gen_sd     = ev_result_per_gen_sd)

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
#' @param pn_gen_sd       genetic standard deviation of trait
#' @param pb_verbose      level of output verbosity
#' @return ev_result      list of two economic values, per unit trait and per genetic sd, if pn_gen_sd is given
compute_ev_price_cont <- function( pn_mean,
                                   pn_sd,
                                   pvec_threshold,
                                   pvec_price,
                                   pn_delta_mean,
                                   pn_gen_sd,
                                   pb_verbose       = FALSE ) {
  ### # getting area under normal distribution for base situation
  vec_freq_base_cum <- pnorm(pvec_threshold, mean = pn_mean, sd = pn_sd)
  vec_freq_base <- diff(c(0,vec_freq_base_cum))
  ### # take last element of vec_freq_base_cum as the cumsum of vec_freq_base
  n_sum_freq_base <- vec_freq_base_cum[length(vec_freq_base_cum)]
  if (n_sum_freq_base > 1){
    stop(" *** Error: compute_ev_price_cont: sum of base frequencies > 1: ", n_sum_freq_base)
  } else {
    vec_freq_base <- c(vec_freq_base, 1-n_sum_freq_base)
  }

  ev_result <- compute_ev_price_disc( pn_mean         = pn_mean,
                                      pn_sd           = pn_sd,
                                      pvec_class_freq = vec_freq_base,
                                      pvec_threshold  = pvec_threshold,
                                      pvec_price      = pvec_price,
                                      pn_delta_mean   = pn_delta_mean,
                                      pn_gen_sd       = pn_gen_sd,
                                      pb_verbose      = pb_verbose)


  return(ev_result)
}


## --- Conversion Functions ----------------------------------------------------
#' @title Relative Economic Factors From Economic Values
#'
#' @description
#' Given a tibble with economic values per genetic standard deviation in
#' parameter ptbl_economic_value, the relative weight of each economic
#' value with respect to the total sum of all economic values is computed.
#'
#' @param ptbl_economic_value tibble with economic values
#' @param pb_first_row_trait_name   flag indicating wether first column are trait names
#' @return tbl_rel_fact tibble with relative factors
#' @export get_relative_economic_factors
get_relative_economic_factors <- function(ptbl_economic_value,
                                          pb_first_row_trait_name = FALSE){

  ### # if pb_first_row_trait_name then remove the first column from ptbl_economic_value
  if (pb_first_row_trait_name){
    mat_economic_value <- as.matrix(ptbl_economic_value[,2:ncol(ptbl_economic_value)])
  } else {
    mat_economic_value <- as.matrix(ptbl_economic_value)
  }

  ### # compute vector of column sums of absolute values of mat_economic_value
  vec_sum_abs_ev <- apply(abs(mat_economic_value), 2, sum)
  ### # extend vector of absolute sums to a matrix
  mat_abs_sum_ev <- matrix(vec_sum_abs_ev, nrow = nrow(mat_economic_value), ncol = ncol(mat_economic_value), byrow = TRUE)
  ### # compute matrix of relative factors
  mat_rel_fact <- mat_economic_value / mat_abs_sum_ev

  ### # check
  # all.equal(sum(apply(abs(mat_factors_ev), 2, sum)),ncol(mat_factors_ev))
  if (! all.equal(sum(apply(abs(mat_rel_fact), 2, sum)), ncol(mat_rel_fact)))
    stop("[ERROR -- get_relative_economic_factors] Sum of relative factors is not 1")

  ### # conversion to tibble
  tbl_rel_fact <- tibble::as_tibble(mat_rel_fact)
  ### # adding first row, if needed
  if (pb_first_row_trait_name){
    tbl_rel_fact <- cbind(ptbl_economic_value[,1], tbl_rel_fact)
  }

  ### # return
  return(tbl_rel_fact)
}
