# METHOD 2



# Betensky R. and Tierney C. (1997). An examination of methods for sample size
# recalculation during an experiment. Statistics in Medicine, 16(22), 2587–2598.
# https://doi.org/10.1002/(SICI)1097-0258(19971130)16:223.0.CO;2-5



# Studying the effect of M value.


# set.seed(048) # Set a seed to ensure reproducible results.



### SECTION NULL: SSR - Betensky, Tierney (where H0 holds) - BLINDED - multiple trial runs ----

# Simulation of a trial under the assumption that H0 holds

 


library(doParallel)
library(foreach)        # Libraries for parallel programming






cores=detectCores()
cores

cl = makeCluster(cores[1]-1)  
registerDoParallel(cl)



M = c(10,20,30,40,50,100,120,150)
rows = length(M)


start_time = Sys.time()


# Executing loop using parallel programming
results = foreach(iter=1:rows, .combine = rbind) %dopar% {                              


  sigmaTRUE = 6.5

  sigma = 6.5 # Assumed parameters

  Delta = 2.34
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 0 # The ACTUAL treatment effect, under the null hypothesis.
  

  n0 = round((((qnorm(1-alpha) + qnorm(1-beta))^2)*(sigma^2)) / (Delta^2)) 
  
  N0 = 2*n0 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N0, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment. 
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  pvalues = rep(NA, Nsim) # An empty vector to store the p-values.
  sigma_reestimates = rep(NA, Nsim) # An empty vector to store the re-estimated sigma values.
  N_reestimates = rep(NA, Nsim) # N re-estimated values
  
  # Loop to perform our simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n0, muC, sigmaTRUE)
    resultsT = rnorm(n0, muT, sigmaTRUE)
    
    interim = round(mean(n0/2))
    
    n_stars = rep(NA,10)
    
    # M times loop.
    for(j in 1:M[iter]){
      
      interim_samples = c(resultsC[1:interim], resultsT[1:interim])  
      sigma_reestimate = sqrt(var(interim_samples))
      
      n = n0
      
      while(n < round((((qnorm(1-alpha)+qnorm(1-beta))^2)*(sigma_reestimate^2)) / (Delta^2)))
      {
        sample_add = sample(interim_samples, size = 1, replace = TRUE)
        
        interim_samples = c(interim_samples, sample_add)
        
        # update n and sigma_reestimate for next iteration of while
        n = length(interim_samples)
        sigma_reestimate = sqrt(var(interim_samples))
        
        cat("\n\nTrial: ",i,"\t\tLoop: ",j)
        cat("\nValue N(both arms): ", 2*n, "\t\tSigma re-estimate: ", sigma_reestimate)
        
      }
      
      n_star = n
      n_stars[j] = n_star
      
    }  
    
    sigma_reestimates[i] = sigma_reestimate  # The final value of sigma at the end of while loop
    
    N_qs = round(mean(n_stars))
    N_reestimates[i] = 2*N_qs
    cat("Trial: ",i,"\t\tN-qsboth arms): ", 2*N_qs)
    
    resultsC_extension = rnorm(N_qs - n0, muC, sigmaTRUE)
    resultsT_extension = rnorm(N_qs - n0, muT, sigmaTRUE)
    
    resultsC_total = c(resultsC, resultsC_extension)
    resultsT_total = c(resultsT, resultsT_extension)
    
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    cat("\nIteration end")
    
    

  }
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  alphahat = mean(rejectnull)
  alphahat
  
  # CI using standard error
  standard_error = sqrt((alphahat*(1 - alphahat))/Nsim)
  standard_error
  
  # Confidence Intervals
  CI = c(alphahat - standard_error, alphahat + standard_error)
  CI
  
  sigma_reestimate = mean(sigma_reestimates)
  sigma_reestimate
  # 500000 simulations: 6.479028
  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))


  
  # N_reestimates          
  
  c(M[iter], sigma_reestimate, round(mean(N_reestimates)), alphahat, CI[1], CI[2])


}


stopCluster(cl)
results

end_time = Sys.time()

met2_h0_null_M = data.frame(results)

names(met2_h0_null_M) = c("M","Sigma Re-estimate", "Sample Size Re-estimated", 
                                "Alpha", "Confidence Interval Low", "Confidence Interval High")

row.names(met2_h0_null_M) = c(1:8)

met2_h0_null_M

write.csv(met2_h0_null_M, "Results/met2_h0_null_M.csv", row.names = FALSE)


print(end_time - start_time)


                                                                             












### SECTION NULL: SSR - Betensky, Tierney (where H1 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H1 holds


cores=detectCores()
cores

cl = makeCluster(cores[1]-1) 
registerDoParallel(cl)



M = c(10,20,30,40,50,100,120,150)
# M = c(10, 100)
rows = length(M)


start_time = Sys.time()

# Executing loop using parallel programming
results2 = foreach(iter=1:rows, .combine = rbind) %dopar% {
  
  
  sigmaTRUE = 6.5
  
  sigma = 6.5 # Assumed parameters
  
  Delta = 2.34 
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 2.34 # The ACTUAL treatment effect, under the alternate hypothesis.
  
  
  # Calculate the desired sample size using the formula (for one group): 
  n0 = round((((qnorm(1-alpha) + qnorm(1-beta))^2)*(sigma^2)) / (Delta^2)) 
  
  N0 = 2*n0 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N0, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment. 
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  pvalues = rep(NA, Nsim) # An empty vector to store the p-values.
  sigma_reestimates = rep(NA, Nsim) # An empty vector to store the re-estimated sigma values.
  N_reestimates = rep(NA, Nsim) # N re-estimated values
  
  # Loop to perform our simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n0, muC, sigmaTRUE)
    resultsT = rnorm(n0, muT, sigmaTRUE)
    
    interim = round(mean(n0/2))
    
    n_stars = rep(NA,10)
    
    # M times loop
    for(j in 1:M[iter]){
      
      interim_samples = c(resultsC[1:interim], resultsT[1:interim])  
      sigma_reestimate = sqrt(var(interim_samples))
      
      n = n0
      
      while(n < round((((qnorm(1-alpha)+qnorm(1-beta))^2)*(sigma_reestimate^2)) / (Delta^2)))
      {
        sample_add = sample(interim_samples, size = 1, replace = TRUE)
        
        interim_samples = c(interim_samples, sample_add)
        
        # update n and sigma_reestimate for next iteration of while
        n = length(interim_samples)
        sigma_reestimate = sqrt(var(interim_samples))
        
        cat("\n\nTrial: ",i,"\t\tLoop: ",j)
        cat("\nValue N(both arms): ", 2*n, "\t\tSigma re-estimate: ", sigma_reestimate)
        
      }
      
      n_star = n
      n_stars[j] = n_star
      
    }  
    
    sigma_reestimates[i] = sigma_reestimate  # The final value of sigma at the end of while loop
    
    N_qs = round(mean(n_stars))
    N_reestimates[i] = 2*N_qs
    cat("Trial: ",i,"\t\tN-qsboth arms): ", 2*N_qs)
    
    resultsC_extension = rnorm(N_qs - n0, muC, sigmaTRUE)
    resultsT_extension = rnorm(N_qs - n0, muT, sigmaTRUE)
    
    resultsC_total = c(resultsC, resultsC_extension)
    resultsT_total = c(resultsT, resultsT_extension)
    
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    cat("\nIteration end")
    
    
    
  }
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  betahat = 1 - mean(rejectnull)
  betahat
  
  # CI using standard error
  standard_error = sqrt((betahat*(1 - betahat))/Nsim)
  standard_error
  
  # Confidence Intervals
  CI = c(betahat - standard_error, betahat + standard_error)
  CI
  
  sigma_reestimate = mean(sigma_reestimates)
  sigma_reestimate
  # 500000 simulations: 6.479028
  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  # N_reestimates
  
  c(M[iter], sigma_reestimate, round(mean(N_reestimates)), betahat, CI[1], CI[2])
  
  
}


stopCluster(cl)
results2

end_time = Sys.time()

met2_h1_null_M = data.frame(results2)

names(met2_h1_null_M) = c("M","Sigma Re-estimate", "Sample Size Re-estimated", 
                          "Beta", "Confidence Interval Beta Low", "Confidence Interval Beta High")

row.names(met2_h1_null_M) = c(1:8)

met2_h1_null_M

write.csv(met2_h1_null_M, "Results/met2_h1_null_M.csv", row.names = FALSE)


print(end_time - start_time)








