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
#' the function compute_asym_selection_response is computing the asymptotic selection response
#'
#' @param
#' @param pvec_economic_values vector of economic values
#' @param pvec_proportion_selected vector indicating the selection intensity
#' @param pvec_generation_interval vector indicating the generation interval
#' @param pvec_number_progeny vector number of offspring
#' @param pvec_proportion_progeny_adult vector number of adult progenies
#' @param pmat_genetic_varcov variance covariance matrix
#' @param pmat_residual_varcov variance covariance residual
#' @param pb_out flag to give out information about the function
#' @return result_asym_sel_res asymptotic selection response
#' @export compute_asym_selection_response
compute_asym_selection_response <- function(pvec_economic_values,
                                            pvec_proportion_selected,
                                            pvec_generation_interval,
                                            pvec_number_progeny,
                                            pvec_proportion_progeny_adult,
                                            pmat_genetic_varcov,
                                            pmat_residual_varcov,
                                            pb_out = FALSE) {
  result_asym_sel_res <- NULL
  n_r_trait <- length(pvec_economic_values)


  ### ##  Proportion of calves
  vec_proportion_calves <- 1 -  pvec_proportion_progeny_adult

  ### ## Starting with male selection candidate
  ### ## Number of offspring per category
  male_number_adults <- floor(pvec_proportion_progeny_adult*pvec_number_progeny[[1]])
  male_number_calves <- pvec_number_progeny[[1]] - male_number_adults



  ### ## Building Design Matrix Z
  ### ## Design Matrix for one calf
  male_designMatrix_calf <- cbind(c(1,0,0),c(0,0,0),c(0,1,0),c(0,0,0),c(0,0,1),c(0,0,0))
  ### ## Design Matrix for one adult
  male_designMatrix_adult <- cbind(c(0,0,0),c(1,0,0),c(0,0,0),c(0,1,0),c(0,0,0),c(0,0,1))
  ### ## Design Matrix for calves
  male_designMatrix_calves <- diag(1,male_number_calves)%x%male_designMatrix_calf
  male_designMatrix_calves <- cbind(male_designMatrix_calves,
                                    matrix(0, nrow=nrow(male_designMatrix_calves),
                                              ncol=male_number_adults*n_r_trait))
  ### ## Design Matrix for adults
  male_designMatrix_adults <- diag(1,male_number_adults)%x%male_designMatrix_adult
  male_designMatrix_adults <- cbind(matrix(0, nrow=nrow(male_designMatrix_adults),
                                           ncol=male_number_calves*n_r_trait),
                                    male_designMatrix_adults)
  ### ## Bind Design Matrix for calves and adults
  male_designMatrixZ <- rbind(male_designMatrix_calves, male_designMatrix_adults)
  ### ## Add selection candidate
  male_designMatrixZ<-cbind(matrix(0, nrow=pvec_number_progeny[[1]]*(n_r_trait/2),ncol=n_r_trait),male_designMatrixZ)
  ### ## Check the dimension of Design Matrix Z
  if(pb_out){
    cat("[compute_asym_selection_response]: check the dimension of Design Matrix Z:\n")
    if(nrow(male_designMatrixZ) != pvec_number_progeny[[1]]*(n_r_trait/2)){
      stop("[compute_asym_selection_response]: something is wrong with the number of rows of Design Matrix Z")
    }
    if(ncol(male_designMatrixZ) != n_r_trait*(pvec_number_progeny[[1]]+1)){
      stop("[compute_asym_selection_response]: something is wrong with the number of columns of Design Matrix Z")
    }
  }



  ### ## Residual Matrix R
  mat_residual_var_cov_calves <- pmat_residual_varcov[c(1,3,5),c(1,3,5)]
  mat_residual_var_cov_adults <- pmat_residual_varcov[c(2,4,6),c(2,4,6)]

  male_calves_kronecker_residual <- diag(male_number_calves) %x% mat_residual_var_cov_calves
  male_adults_kronecker_residual <- diag(male_number_adults) %x% mat_residual_var_cov_adults
  male_calves_kronecker_residual_extended <- cbind(male_calves_kronecker_residual,
                                                   matrix(0,nrow=nrow(male_calves_kronecker_residual),
                                                            ncol=ncol(male_adults_kronecker_residual)))
  male_adults_kronecker_residual_extended <- cbind(matrix(0,nrow=nrow(male_adults_kronecker_residual),
                                                             ncol=ncol(male_calves_kronecker_residual)),
                                                  male_adults_kronecker_residual)
  male_ResidualMatrix <- rbind(male_calves_kronecker_residual_extended,male_adults_kronecker_residual_extended)




  ### ## Relationship Matrix A
  mnumb <- pvec_number_progeny[[1]]+1
  male_Pedigree <- pedigreemm::pedigree(sire = c(NA,rep(1,times = pvec_number_progeny[[1]])),
                                        dam =c(NA,rep(NA,times = pvec_number_progeny[[1]])),
                                        label = 1:mnumb)
  ### ## Compute inverse of A
  male_AInv <- as.matrix(pedigreemm::getAInv(male_Pedigree))



  ### ## Computation of PEV
  male_PEV <- solve(t(male_designMatrixZ)%*%solve(male_ResidualMatrix)%*%male_designMatrixZ+male_AInv%x%solve(pmat_genetic_varcov))
  if(pb_out){
    cat("[compute_asym_selection_response]: check the dimension of PEV Z:\n")
    if(nrow(male_PEV) != n_r_trait*(pvec_number_progeny[[1]]+1)){
      stop("[compute_asym_selection_response]: something is wrong with the number of rows of PEV")
    }
    if(ncol(male_PEV) != n_r_trait*(pvec_number_progeny[[1]]+1)){
      stop("[compute_asym_selection_response]: something is wrong with the number of columns of PEV")
    }
  }
  ### ## Computation of PEV for selection candidate
  male_PEV <- male_PEV[1:n_r_trait,1:n_r_trait]



  ### ## Computation of variance covariance matrix for predicting breeding values
  mat_male_var_cov_pbv <- pmat_genetic_varcov - male_PEV




  ### ## -----------------------------------------
  ### ## Following with female selection candidate
  ### ## Number of offspring per category
  female_number_adults <- floor(pvec_proportion_progeny_adult*pvec_number_progeny[[2]])
  female_number_calves <- pvec_number_progeny[[2]] - female_number_adults



  ### ## Building Design Matrix Z
  ### ## Design Matrix for one calf
  female_designMatrix_calf <- cbind(c(1,0,0),c(0,0,0),c(0,1,0),c(0,0,0),c(0,0,1),c(0,0,0))
  ### ## Design Matrix for one adult
  female_designMatrix_adult <- cbind(c(0,0,0),c(1,0,0),c(0,0,0),c(0,1,0),c(0,0,0),c(0,0,1))
  ### ## Design Matrix for calves
  female_designMatrix_calves <- diag(1,female_number_calves)%x%female_designMatrix_calf
  female_designMatrix_calves <- cbind(female_designMatrix_calves,
                                    matrix(0, nrow=nrow(female_designMatrix_calves),
                                           ncol=female_number_adults*n_r_trait))
  ### ## Design Matrix for adults
  female_designMatrix_adults <- diag(1,female_number_adults)%x%female_designMatrix_adult
  female_designMatrix_adults <- cbind(matrix(0, nrow=nrow(female_designMatrix_adults),
                                           ncol=female_number_calves*n_r_trait),
                                      female_designMatrix_adults)
  ### ## Bind Design Matrix for calves and adults
  female_designMatrixZ <- rbind(female_designMatrix_calves, female_designMatrix_adults)
  ### ## Add selection candidate
  female_designMatrixZ<-cbind(matrix(0, nrow=pvec_number_progeny[[2]]*(n_r_trait/2),ncol=n_r_trait),female_designMatrixZ)
  ### ## Check the dimension of Design Matrix Z
  if(pb_out){
    cat("[compute_asym_selection_response]: check the dimension of Design Matrix Z:\n")
    if(nrow(female_designMatrixZ) != pvec_number_progeny[[2]]*(n_r_trait/2)){
      stop("[compute_asym_selection_response]: something is wrong with the number of rows of Design Matrix Z")
    }
    if(ncol(female_designMatrixZ) != n_r_trait*(pvec_number_progeny[[2]]+1)){
      stop("[compute_asym_selection_response]: something is wrong with the number of columns of Design Matrix Z")
    }
  }



  ### ## Residual Matrix R
  mat_residual_var_cov_calves <- pmat_residual_varcov[c(1,3,5),c(1,3,5)]
  mat_residual_var_cov_adults <- pmat_residual_varcov[c(2,4,6),c(2,4,6)]

  female_calves_kronecker_residual <- diag(female_number_calves) %x% mat_residual_var_cov_calves
  female_adults_kronecker_residual <- diag(female_number_adults) %x% mat_residual_var_cov_adults
  female_calves_kronecker_residual_extended <- cbind(female_calves_kronecker_residual,
                                                   matrix(0,nrow=nrow(female_calves_kronecker_residual),
                                                          ncol=ncol(female_adults_kronecker_residual)))
  female_adults_kronecker_residual_extended <- cbind(matrix(0,nrow=nrow(female_adults_kronecker_residual),
                                                          ncol=ncol(female_calves_kronecker_residual)),
                                                   female_adults_kronecker_residual)
  female_ResidualMatrix <- rbind(female_calves_kronecker_residual_extended,female_adults_kronecker_residual_extended)




  ### ## Relationship Matrix A
  fnumb <- pvec_number_progeny[[2]]+1
  female_Pedigree <- pedigreemm::pedigree(sire = c(NA,rep(NA,times = pvec_number_progeny[[2]])),
                                          dam = c(NA,rep(1,times = pvec_number_progeny[[2]])),
                                         label = 1:fnumb)
  ### ## Compute inverse of A
  female_AInv <- as.matrix(pedigreemm::getAInv(female_Pedigree))



  ### ## Computation of PEV
  female_PEV <- solve(t(female_designMatrixZ)%*%solve(female_ResidualMatrix)%*%female_designMatrixZ+female_AInv%x%solve(pmat_genetic_varcov))
  if(pb_out){
    cat("[compute_asym_selection_response]: check the dimension of PEV Z:\n")
    if(nrow(female_PEV) != n_r_trait*(pvec_number_progeny[[2]]+1)){
      stop("[compute_asym_selection_response]: something is wrong with the number of rows of PEV")
    }
    if(ncol(female_PEV) != n_r_trait*(pvec_number_progeny[[2]]+1)){
      stop("[compute_asym_selection_response]: something is wrong with the number of columns of PEV")
    }
  }
  ### ## Computation of PEV for selection candidate
  female_PEV <- female_PEV[1:n_r_trait,1:n_r_trait]



  ### ## Computation of variance covariance matrix for predicting breeding values
  mat_female_var_cov_pbv <- pmat_genetic_varcov - female_PEV





  ### ## -----------------------------------------
  ### ## Selection intensity and Factor of variance reduction k
  male_i <- dnorm(qnorm(1-pvec_proportion_selected[[1]]))/pvec_proportion_selected[[1]]
  male_x <- qnorm(pvec_proportion_selected[[1]], lower.tail = FALSE)
  male_k <- male_i*(male_i-male_x)

  female_i <- dnorm(qnorm(1-pvec_proportion_selected[[2]]))/pvec_proportion_selected[[2]]
  female_x <- qnorm(pvec_proportion_selected[[2]], lower.tail = FALSE)
  female_k <- female_i*(female_i-female_x)

  ### ## Computation of asymptotic variance covariance matrix
  mat_kInv <- solve(matrix(c(1+0.5*male_k, 0.5*female_k,0.5*male_k, 1+0.5*female_k),nrow = 2, ncol = 2, byrow=TRUE))
  mat_var_fact_reduction <- mat_kInv%x%diag(1,n_r_trait)

  mat_var_cov_pbv <- rbind(mat_male_var_cov_pbv, mat_female_var_cov_pbv)
  mat_asympt_var_cov_pbv <- mat_var_fact_reduction%*%mat_var_cov_pbv





  ### ## -----------------------------------------
  ### ## Asymptotic genetic gain of the estimated overall breeding value (SZENARIO 2)
  ### ## Separate asymptotic variance covariance matrix for male and female
  mat_male_asympt_var_cov_pbv <- mat_asympt_var_cov_pbv[1:n_r_trait,]
  mat_female_asympt_var_cov_pbv <- mat_asympt_var_cov_pbv[(n_r_trait+1):nrow(mat_asympt_var_cov_pbv),]

  result_asym_sel_res <- (male_i*sqrt(t(as.matrix(pvec_economic_values))%*%mat_male_asympt_var_cov_pbv%*%as.matrix(pvec_economic_values))+
                          female_i*sqrt(t(as.matrix(pvec_economic_values))%*%mat_female_asympt_var_cov_pbv%*%as.matrix(pvec_economic_values)))/
                          sum(pvec_generation_interval)


  return(result_asym_sel_res)
}
