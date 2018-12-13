###
###
###
###   Purpose: Functions To Compute Asymptotic Selection Response
###   started: 2018-12-13 (skn and pvr)
###
### ################################################################ ###


#' @title Compute asymptotic selection response
#'
#' @description
#'
#' @param
#'
compute_asym_selection_response <- function(pvec_economic_values,
                                            pvec_proportion_selected,
                                            pvec_generation_interval,
                                            pvec_number_progeny,
                                            pvec_proportion_progeny_adult,
                                            pmat_genetic_varcov,
                                            pmat_residual_varcov) {
  result_asym_sel_res <- NULL

  ### ##  Proportion of calves
  vec_proportion_calves <- 1 -  pvec_proportion_progeny_adult

  ### ## Starting with male selection candidate
  ### ## Number of offspring per category
  male_number_adults <- floor(pvec_proportion_progeny_adult[[1]]*pvec_number_progeny[[1]])
  male_number_calves <- pvec_number_progeny[[1]] - male_number_adults

  ### ##










  return(result_asym_sel_res)
}
