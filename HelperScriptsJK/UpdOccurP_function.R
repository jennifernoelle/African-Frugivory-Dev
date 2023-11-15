#' Function that updates the probability of occurence. 
#' Previously occurrence probabilities were specified via the arguments occur_B, occur_P
#' Now prior probabilities will be optionally supplied via the arguments p_occur_B, p_occur_P
#' We also add a new argument to the sampling vector to indicate whether or not these should be sampled
#' If not: we assume occurrence probabilities are known
#' If so: we assume occurrence probablities are known with uncertainty
#' 
#' 
#' Each species might occur in the study area or not. Even though probabilities
#' of occurence can be provided, the detected interactions also inform us of
#' whether the species occur. QUESTION: does this also hold for occurrence 
#' probabilities?
#' 
#' @param detected Matrix. Rows are species, columns are studies. Indicator of
#' whether the species were detected to be present in the given study. Derived
#' from the provided occurrence matrix, and checked against obs_A. Note that
#' while most studies only report interacting species, some studies do report
#' present but not interacting species. This often occurs when subsetting. 
#' @param occur Matrix. Rows are species, columns are studies. It includes the
#' prior probabilities of occurrence.
#' @param occor_others. Matrix Rows are species, columns are studies. It includes
#' the current values of occurrence for the other set of species, the ones that
#' are not being updated. 
#' @param probobs_curr Vector of current values for the probability of
#' observation for the set of species we are updating.
#' @param probobs_others Vector of current values for the probability of
#' observation for the other set of species, the ones that are not being
#' updated.
#' @param curr_inter Matrix with 0,1 entries for the current posterior samples
#' of possible interactions. The rows correspond to the species we are updating
#' and the columns correspond to the other set of species.
#' @param focus Array of three dimensions corresponding to the two sets of
#' species, and the different studies. Values are 0 or 1. This array represents
#' whether the current study would be willing to record the interactions of a
#' given species. The value should be 1 except for studies that are animal or
#' plant-oriented, where the value for species that are not of focus will be 0.
#' 
#UpdOccurP <- function(detected, occur, occur_others, probobs_curr, probobs_others, curr_inter,
#                     focus) {
  
# Args it needs: 
# detection indicator: detected_B, detected_P
# current value of occurrence indicator for this species type only: this_O_B, this_O_P
# Prior prob of occurrence: come up with varname: prior_occur_B, prior_occur_P (prior mean)
# Prior mean for occurrence: default is 1, varname: mh_occ_sd
# current value of occurrence probability: occur_B, occur_P (old version: these are fixed, now we're updating)
# MH step size: mh_occ_step

# Other function things to do: add argument, create default assignment for prior probs and default mh step size and prior sd

#' Function call will be like this

# # Line by line debugging
# occur_indicators = this_O_B
# occur_prior_probs = prior_occur_B = occur_B # just for now because I haven't set the new argument yet
# curr_occur_probs = occur_B
# mh_occ_sd = 0.1
# mh_occ_step = 0.25
# p_occ_accepted = matrix(0, nrow = nrow(this_O_B), ncol = ncol(this_O_B)) 
# detected <- detected_B
# 
# 
# upd_P_O_B <- UpdOccurP(occur_indicators = this_O_B, 
#                         occur_prior_probs = prior_occur_B, 
#                         curr_occur_probs = occur_B, 
#                         detected = detected_B,
#                         mh_occ_sd = mh_occ_sd, 
#                         mh_occ_step = mh_occ_step)
# # Then in the mcmc function next step is to update the proobs
# this_occur_B <- upd_P_O_B$new_value_p


UpdOccurP <- function(occur_indicators, occur_prior_probs, curr_occur_probs, detected, mh_occ_sd, mh_occ_step) {
  
  num_obs <- nrow(occur_indicators)
  num_studies <- ncol(occur_indicators)
  
  # Return something like this to keep track of acceptance?
  new_value_p = curr_occur_probs
  accepted = matrix(0, nrow = num_obs, ncol = num_studies)
  
  # If the species were observed to interact in the study, the occurrence probability
  # is always equal to 1. First, I will draw the occurrence probabiilty
  # as if the species were not observed to interact. Then, I will assign the
  # ones observed interacting to 1.
  
  # Updating the occurrence probability as if the species was not detected in
  # the study. For all species, one study at a time.
  
  
  for (st in 1 : num_studies) {
    
    # Generating the proposal
    occur <- occur_indicators[, st]
    p_curr <- curr_occur_probs[, st]
    p_prior <- occur_prior_probs[, st]
    #p_prop <- rtnorm(n = num_obs,  mean = p_curr, sd = mh_occ_step, lower = 0, upper = 1)
    p_prop <- rtruncnorm(n = 1,  mean = p_curr, sd = mh_occ_step, a = 0, b = 1)

    # Calculating the likelihood, prior, and proposal for proposed value for each proposed value: 
    #log_prior_prop <- dtnorm(x = p_prop, lower = 0, upper = 1, mean = p_prior, sd = mh_occ_sd, log = TRUE)
    log_prior_prop <- log(dtruncnorm(x = p_prop, a = 0, b = 1, mean = p_prior, sd = mh_occ_sd))
    log_lik_prop <- log(p_prop^occur * (1-p_prop)^(1-occur))
    #log_prop <- dtnorm(x = p_prop, lower = 0, upper = 1, mean = p_curr, sd = mh_occ_sd, log = TRUE) # msm version is not significantly faster
    log_prop <- log(dtruncnorm(x = p_prop, a = 0, b= 1, mean = p_curr, sd = mh_occ_sd))  
    
    # Calculating the likelihood, prior, and proposal for proposed value for each current value: 
    #log_prior_curr <- dtnorm(x = p_curr, lower = 0, upper = 1, mean = p_prior, sd = mh_occ_sd, log = TRUE)
    log_prior_curr <- log(dtruncnorm(x = p_curr, a = 0, b = 1, mean = p_prior, sd = mh_occ_sd))
    log_lik_curr <- log(p_curr^occur * (1-p_curr)^(1-occur))
    #log_curr <- dtnorm(x = p_curr, lower = 0, upper = 1, mean = p_prop, sd = mh_occ_sd, log = TRUE) 
    log_curr <- log(dtruncnorm(x = p_curr, a = 0, b = 1, mean = p_prop, sd = mh_occ_sd))
    
    AP <- log_lik_prop + log_prior_prop + log_curr
    AP <- AP - (log_lik_curr + log_prior_curr + log_prop)
    
    # 
    # AP_prop <- log_lik_prop + log_prior_prop + log_curr
    # AP_curr <- log_lik_curr + log_prior_curr + log_prop
    # 
    update <- log(runif(num_obs)) < AP
    new_value_p[update, st] <- p_prop[update]
    accepted[update, st] <- 1

  }
  
  # Now, if a species was detected to interact in a study, set their occurrence
  # probability to 1.
  wh_detected <- which(detected == 1)
  new_value_p[wh_detected] <- 1
  
  r <- list(new_value_p = new_value_p, accepted = accepted)
  return(r)
  
  # Debugging
  # not_detected <- which(detected == 0)
  # hist(new_value_p[not_detected])
}

## debug: why is acceptance rate so high? Maybe because i'm starting at prior mean?
## debug: it's actually possible for a plant to be recorded as present without detected interactions
# use prior == 1 instead of detected ==1 

