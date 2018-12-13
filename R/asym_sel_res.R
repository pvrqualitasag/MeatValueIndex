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
                                            pmat_residual_varcov,
                                            pb_out = FALSE) {
  result_asym_sel_res <- NULL
  n_r_trait <- length(pvec_economic_values)


  ### ##  Proportion of calves
  vec_proportion_calves <- 1 -  pvec_proportion_progeny_adult

  ### ## Starting with male selection candidate
  ### ## Number of offspring per category
  male_number_adults <- floor(pvec_proportion_progeny_adult[[1]]*pvec_number_progeny[[1]])
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
  male_designMatrix_adults <- cbind(matrix(0, nrow=nrow(male_designMatrix_adult),
                                           ncol=male_number_calves*n_r_trait),
                                    male_designMatrix_adult)
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
  male_PEV <- male_PEV[1:n,1:n]



  ### ## Computation of variance covariance matrix for predicting breeding values
  mat_male_var_cov_pbv <- pmat_genetic_varcov - male_PEV










  return(result_asym_sel_res)
}
