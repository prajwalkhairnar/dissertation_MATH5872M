# Studying the effect of change in allocation ratio (k)


# set.seed(048) # Setting a seed to ensure reproducible results.


## Null Hypothesis H0 Holds  ----

### SECTION I: SSR - Chow, Wang and Shao - k = 0.25 - (where H0 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H0 holds


iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )   


rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_alpha_reestimate = rep(NA, rows)
iter_CI_low = rep(NA, rows)
iter_CI_high = rep(NA, rows)

df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_alpha_reestimate = rep(NA, rows)
non_ssr_iter_CI_low  = rep(NA, rows)
non_ssr_iter_CI_high = rep(NA, rows)

for(iter in 1:rows)
{
  
  sigmaTRUE = iter_sigmaTRUE[iter]
  
  sigma = iter_sigma[iter]

  Delta = 2.34 
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 0 # The ACTUAL treatment effect, under the null hypothesis.
  
  
  # Calculate the desired sample size using the formula (for one group) using n2: 
  
  k = 0.25  # Allocation ratio 
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   (based on delta = 0, should be same. NO DIFFERENCE)
  
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  
  pvalues = rep(NA, Nsim) # An empty vector to store the p-values.
  sigma_reestimates = rep(NA, Nsim) # An empty vector to store the re-estimated sigma values.
  N_reestimates = rep(NA, Nsim) # N re-estimated values
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  # Loop to perform simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])  
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                    
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    # t-test on re-estimated sample size data
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
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

  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_alpha_reestimate[iter] = alphahat
  
  iter_CI_low[iter] = CI[1]
  iter_CI_high[iter] = CI[2]
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_alphahat = mean(non_ssr_rejectnull)
  
  non_ssr_standard_error = sqrt((non_ssr_alphahat*(1 - non_ssr_alphahat))/Nsim)
  non_ssr_standard_error
  
  non_ssr_CI = c(non_ssr_alphahat - non_ssr_standard_error, non_ssr_alphahat + non_ssr_standard_error)
  
  non_ssr_iter_alpha_reestimate[iter] = non_ssr_alphahat
  
  non_ssr_iter_CI_low[iter] = non_ssr_CI[1]
  non_ssr_iter_CI_high[iter] = non_ssr_CI[2]
  
  
}



results_k0.25 = data.frame("delta" = rep(0,rows),
                           "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                           "Sigma Re-estimate" = iter_sigma_reestimate, 
                           "Sample Size Initial" = non_ssr_ss , "Sample Size Re-estimate" = iter_ss_reestimate, 
                           "Alpha Re-estimate" = iter_alpha_reestimate,
                           "Confidence Interval Low" = iter_CI_low, "Confidence Interval High" = iter_CI_high, 
                           "NON SSR Alpha Re-estimate" = non_ssr_iter_alpha_reestimate, 
                           "NON SSR Confidence Interval Low" = non_ssr_iter_CI_low, "NON SSR Confidence Interval High" = non_ssr_iter_CI_high)
print(results_k0.25)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))


# Confidence Intervals for the re-estimated sample size
ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_k0.25$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_k0.25$Sample.Size.Re.estimate = ss_reestimate


write.csv(results_k0.25, "Results/allocation_study_blinded_h0_k_0;25.csv", row.names = FALSE)


write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h0_k_0;25_distribution.csv", row.names = FALSE)

































### SECTION II: SSR - Chow, Wang and Shao - k = 0.50 - (where H0 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H0 holds


iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )  

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_alpha_reestimate = rep(NA, rows)
iter_CI_low = rep(NA, rows)
iter_CI_high = rep(NA, rows)

df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_alpha_reestimate = rep(NA, rows)
non_ssr_iter_CI_low  = rep(NA, rows)
non_ssr_iter_CI_high = rep(NA, rows)

for(iter in 1:rows)
{

  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters
  sigma = iter_sigma[iter]
  
  
  
  Delta = 2.34 
  
  alpha = 0.05 # Desired error rates for 
  beta = 0.2   # simulated trial.
  
  delta = 0 # The ACTUAL treatment effect, under the null hypothesis.
  
  

  k = 0.5  # Allocation ratio 
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   (based on delta = 0, should be same. NO DIFFERENCE)
  
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) # An empty vector to store the p-values.
  sigma_reestimates = rep(NA, Nsim) # An empty vector to store the re-estimated sigma values.
  N_reestimates = rep(NA, Nsim) # N re-estimated values
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  # Loop to perform our simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
   
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])  
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                    
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
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

  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_alpha_reestimate[iter] = alphahat
  
  iter_CI_low[iter] = CI[1]
  iter_CI_high[iter] = CI[2]
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_alphahat = mean(non_ssr_rejectnull)
  
  non_ssr_standard_error = sqrt((non_ssr_alphahat*(1 - non_ssr_alphahat))/Nsim)
  non_ssr_standard_error
  
  non_ssr_CI = c(non_ssr_alphahat - non_ssr_standard_error, non_ssr_alphahat + non_ssr_standard_error)
  
  non_ssr_iter_alpha_reestimate[iter] = non_ssr_alphahat
  
  non_ssr_iter_CI_low[iter] = non_ssr_CI[1]
  non_ssr_iter_CI_high[iter] = non_ssr_CI[2]
  
  
}



results_k0.50 = data.frame("delta" = rep(0,rows),
                           "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                           "Sigma Re-estimate" = iter_sigma_reestimate, 
                           "Sample Size Initial" = non_ssr_ss , "Sample Size Re-estimate" = iter_ss_reestimate, 
                           "Alpha Re-estimate" = iter_alpha_reestimate,
                           "Confidence Interval Low" = iter_CI_low, "Confidence Interval High" = iter_CI_high, 
                           "NON SSR Alpha Re-estimate" = non_ssr_iter_alpha_reestimate, 
                           "NON SSR Confidence Interval Low" = non_ssr_iter_CI_low, "NON SSR Confidence Interval High" = non_ssr_iter_CI_high)
print(results_k0.50)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_k0.50$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_k0.50$Sample.Size.Re.estimate = ss_reestimate

write.csv(results_k0.50, "Results/allocation_study_blinded_h0_k_0;50.csv", row.names = FALSE)


write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h0_k_0;50_distribution.csv", row.names = FALSE)
































### SECTION III: SSR - Chow, Wang and Shao - k = 0.75 - (where H0 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H0 holds



iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )   

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_alpha_reestimate = rep(NA, rows)
iter_CI_low = rep(NA, rows)
iter_CI_high = rep(NA, rows)


df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_alpha_reestimate = rep(NA, rows)
non_ssr_iter_CI_low  = rep(NA, rows)
non_ssr_iter_CI_high = rep(NA, rows)

for(iter in 1:rows)
{
  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters
  sigma = iter_sigma[iter]
  
  
  
  Delta = 2.34 
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 0 # The ACTUAL treatment effect, under the null hypothesis.
  
  

  k = 0.75  # Allocation ratio
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   (based on delta = 0, should be same. NO DIFFERENCE)
  
  
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) # An empty vector to store the p-values.
  sigma_reestimates = rep(NA, Nsim) # An empty vector to store the re-estimated sigma values.
  N_reestimates = rep(NA, Nsim) # N re-estimated values
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  # Loop to perform our simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])  
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                    
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                    
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
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
````
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_alpha_reestimate[iter] = alphahat
  
  iter_CI_low[iter] = CI[1]
  iter_CI_high[iter] = CI[2]
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_alphahat = mean(non_ssr_rejectnull)
  
  non_ssr_standard_error = sqrt((non_ssr_alphahat*(1 - non_ssr_alphahat))/Nsim)
  non_ssr_standard_error
  
  non_ssr_CI = c(non_ssr_alphahat - non_ssr_standard_error, non_ssr_alphahat + non_ssr_standard_error)
  
  non_ssr_iter_alpha_reestimate[iter] = non_ssr_alphahat
  
  non_ssr_iter_CI_low[iter] = non_ssr_CI[1]
  non_ssr_iter_CI_high[iter] = non_ssr_CI[2]
  
}



results_k0.75 = data.frame("delta" = rep(0,rows),
                           "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                           "Sigma Re-estimate" = iter_sigma_reestimate, 
                           "Sample Size Initial" = non_ssr_ss , "Sample Size Re-estimate" = iter_ss_reestimate, 
                           "Alpha Re-estimate" = iter_alpha_reestimate,
                           "Confidence Interval Low" = iter_CI_low, "Confidence Interval High" = iter_CI_high, 
                           "NON SSR Alpha Re-estimate" = non_ssr_iter_alpha_reestimate, 
                           "NON SSR Confidence Interval Low" = non_ssr_iter_CI_low, "NON SSR Confidence Interval High" = non_ssr_iter_CI_high)

print(results_k0.75)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_k0.75$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_k0.75$Sample.Size.Re.estimate = ss_reestimate

write.csv(results_k0.75, "Results/allocation_study_blinded_h0_k_0;75.csv", row.names = FALSE)


write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h0_k_0;75_distribution.csv", row.names = FALSE)


































### SECTION IV: SSR - Chow, Wang and Shao - k = 1 - (where H0 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H0 holds



iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )  

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_alpha_reestimate = rep(NA, rows)
iter_CI_low = rep(NA, rows)
iter_CI_high = rep(NA, rows)

df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_alpha_reestimate = rep(NA, rows)
non_ssr_iter_CI_low  = rep(NA, rows)
non_ssr_iter_CI_high = rep(NA, rows)

for(iter in 1:rows)
{

  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters for us to plug
  sigma = iter_sigma[iter]
  
  
  
  Delta = 2.34
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 0 # The ACTUAL treatment effect, under the null hypothesis.
  
  

  k = 1  # Allocation ratio
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   (based on delta = 0, should be same. NO DIFFERENCE)
  
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim)
  sigma_reestimates = rep(NA, Nsim) 
  N_reestimates = rep(NA, Nsim) 
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
 
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])  
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                     
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    # t-test on re-estimated sample size data
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
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

    
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_alpha_reestimate[iter] = alphahat
  
  iter_CI_low[iter] = CI[1]
  iter_CI_high[iter] = CI[2]
  
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_alphahat = mean(non_ssr_rejectnull)
  
  non_ssr_standard_error = sqrt((non_ssr_alphahat*(1 - non_ssr_alphahat))/Nsim)
  non_ssr_standard_error
  
  non_ssr_CI = c(non_ssr_alphahat - non_ssr_standard_error, non_ssr_alphahat + non_ssr_standard_error)
  
  non_ssr_iter_alpha_reestimate[iter] = non_ssr_alphahat
  
  non_ssr_iter_CI_low[iter] = non_ssr_CI[1]
  non_ssr_iter_CI_high[iter] = non_ssr_CI[2]
  
  
}



results_k1.00 = data.frame("delta" = rep(0,rows),
                           "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                           "Sigma Re-estimate" = iter_sigma_reestimate, 
                           "Sample Size Initial" = non_ssr_ss , "Sample Size Re-estimate" = iter_ss_reestimate, 
                           "Alpha Re-estimate" = iter_alpha_reestimate,
                           "Confidence Interval Low" = iter_CI_low, "Confidence Interval High" = iter_CI_high, 
                           "NON SSR Alpha Re-estimate" = non_ssr_iter_alpha_reestimate, 
                           "NON SSR Confidence Interval Low" = non_ssr_iter_CI_low, "NON SSR Confidence Interval High" = non_ssr_iter_CI_high)

print(results_k1.00)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_k1.00$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_k1.00$Sample.Size.Re.estimate = ss_reestimate

write.csv(results_k1.00, "Results/allocation_study_blinded_h0_k_1;00.csv", row.names = FALSE)


write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h0_k_1;00_distribution.csv", row.names = FALSE)







































### SECTION V: SSR - Chow, Wang and Shao - k = 1.5 - (where H0 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H0 holds



iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )   

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_alpha_reestimate = rep(NA, rows)
iter_CI_low = rep(NA, rows)
iter_CI_high = rep(NA, rows)



df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_alpha_reestimate = rep(NA, rows)
non_ssr_iter_CI_low  = rep(NA, rows)
non_ssr_iter_CI_high = rep(NA, rows)



for(iter in 1:rows)
{

  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters for us to plug
  sigma = iter_sigma[iter]
  
  
  
  Delta = 2.34 
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 0 # The ACTUAL treatment effect, under the null hypothesis.
  
  

  k = 1.5  # Allocation ratio 
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   (based on delta = 0, should be same. NO DIFFERENCE)
  
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) 
  sigma_reestimates = rep(NA, Nsim) 
  N_reestimates = rep(NA, Nsim)
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
 
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])  
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                     
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                    
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    # t-test on re-estimated sample size data
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
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

  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_alpha_reestimate[iter] = alphahat
  
  iter_CI_low[iter] = CI[1]
  iter_CI_high[iter] = CI[2]
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_alphahat = mean(non_ssr_rejectnull)
  
  non_ssr_standard_error = sqrt((non_ssr_alphahat*(1 - non_ssr_alphahat))/Nsim)
  non_ssr_standard_error
  
  non_ssr_CI = c(non_ssr_alphahat - non_ssr_standard_error, non_ssr_alphahat + non_ssr_standard_error)
  
  non_ssr_iter_alpha_reestimate[iter] = non_ssr_alphahat
  
  non_ssr_iter_CI_low[iter] = non_ssr_CI[1]
  non_ssr_iter_CI_high[iter] = non_ssr_CI[2]
  
}



results_k1.50 = data.frame("delta" = rep(0,rows),
                           "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                           "Sigma Re-estimate" = iter_sigma_reestimate, 
                           "Sample Size Initial" = non_ssr_ss , "Sample Size Re-estimate" = iter_ss_reestimate, 
                           "Alpha Re-estimate" = iter_alpha_reestimate,
                           "Confidence Interval Low" = iter_CI_low, "Confidence Interval High" = iter_CI_high, 
                           "NON SSR Alpha Re-estimate" = non_ssr_iter_alpha_reestimate, 
                           "NON SSR Confidence Interval Low" = non_ssr_iter_CI_low, "NON SSR Confidence Interval High" = non_ssr_iter_CI_high)

print(results_k1.50)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_k1.50$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_k1.50$Sample.Size.Re.estimate = ss_reestimate

write.csv(results_k1.50, "Results/allocation_study_blinded_h0_k_1;50.csv", row.names = FALSE)


write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h0_k_1;50_distribution.csv", row.names = FALSE)


































### SECTION VI: SSR - Chow, Wang and Shao - k = 2 - (where H0 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H0 holds


iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 ) 

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_alpha_reestimate = rep(NA, rows)
iter_CI_low = rep(NA, rows)
iter_CI_high = rep(NA, rows)


df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_alpha_reestimate = rep(NA, rows)
non_ssr_iter_CI_low  = rep(NA, rows)
non_ssr_iter_CI_high = rep(NA, rows)

for(iter in 1:rows)
{

  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters for us to plug
  sigma = iter_sigma[iter]
  
  
  
  Delta = 2.34 
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 0 # The ACTUAL treatment effect, under the null hypothesis.
  
  
  k = 2  # Allocation ratio
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   (based on delta = 0, should be same. NO DIFFERENCE)
  
  
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) 
  sigma_reestimates = rep(NA, Nsim)
  N_reestimates = rep(NA, Nsim)
  
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  # Loop to perform our simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
  
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])  
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                     
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    # t-test on re-estimated sample size data
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
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

  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_alpha_reestimate[iter] = alphahat
  
  iter_CI_low[iter] = CI[1]
  iter_CI_high[iter] = CI[2]
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_alphahat = mean(non_ssr_rejectnull)
  
  non_ssr_standard_error = sqrt((non_ssr_alphahat*(1 - non_ssr_alphahat))/Nsim)
  non_ssr_standard_error
  
  non_ssr_CI = c(non_ssr_alphahat - non_ssr_standard_error, non_ssr_alphahat + non_ssr_standard_error)
  
  non_ssr_iter_alpha_reestimate[iter] = non_ssr_alphahat
  
  non_ssr_iter_CI_low[iter] = non_ssr_CI[1]
  non_ssr_iter_CI_high[iter] = non_ssr_CI[2]
  
}



results_k2.00 = data.frame("delta" = rep(0,rows),
                           "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                           "Sigma Re-estimate" = iter_sigma_reestimate, 
                           "Sample Size Initial" = non_ssr_ss , "Sample Size Re-estimate" = iter_ss_reestimate, 
                           "Alpha Re-estimate" = iter_alpha_reestimate,
                           "Confidence Interval Low" = iter_CI_low, "Confidence Interval High" = iter_CI_high, 
                           "NON SSR Alpha Re-estimate" = non_ssr_iter_alpha_reestimate, 
                           "NON SSR Confidence Interval Low" = non_ssr_iter_CI_low, "NON SSR Confidence Interval High" = non_ssr_iter_CI_high)

print(results_k2.00)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_k2.00$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_k2.00$Sample.Size.Re.estimate = ss_reestimate

write.csv(results_k2.00, "Results/allocation_study_blinded_h0_k_2;00.csv", row.names = FALSE)


write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h0_k_2;00_distribution.csv", row.names = FALSE)





































## Alternate Hypothesis H1 Holds ----



### SECTION I: SSR - Chow, Wang and Shao - k = 0.25 - (where H1 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H1 holds







iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 ) 

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_beta_reestimate = rep(NA, rows)
iter_power_reestimate = rep(NA, rows)


iter_CI_beta_low = rep(NA, rows)
iter_CI_beta_high = rep(NA, rows)
iter_CI_power_low = rep(NA, rows)
iter_CI_power_high = rep(NA, rows)



df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_beta_reestimate = rep(NA, rows)
non_ssr_iter_power_reestimate = rep(NA, rows)

non_ssr_iter_CI_beta_low  = rep(NA, rows)
non_ssr_iter_CI_beta_high = rep(NA, rows)

non_ssr_iter_CI_power_low  = rep(NA, rows)
non_ssr_iter_CI_power_high = rep(NA, rows)

for(iter in 1:rows)
{

  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters for us to plug
  sigma = iter_sigma[iter]
  
  
  Delta = 2.34
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 2.34 # The ACTUAL treatment effect, under the alternate hypothesis.
  
  

  k = 0.25  # Allocation ratio 
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.
  
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim)
  sigma_reestimates = rep(NA, Nsim) 
  N_reestimates = rep(NA, Nsim) 
  
  non_ssr_pvalues = rep(NA, Nsim)
  

  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
   
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    
    
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                     
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    # t-test on re-estimated sample size data
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
  betahat = 1 - mean(rejectnull)
  betahat

  standard_error_1 = sqrt((betahat*(1 - betahat))/Nsim)
  standard_error_1
  
  
  # Confidence Intervals
  CI_beta = c(betahat - standard_error_1, betahat + standard_error_1)
  CI_beta
  
  
  sigma_reestimate = mean(sigma_reestimates)
  sigma_reestimate

  power = 1 - betahat
  cat("Power : ",power)
  
  standard_error_2 = sqrt((power*(1 - power))/Nsim)
  standard_error_2
  
  
  # Confidence Intervals
  CI_power = c(power - standard_error_2, power + standard_error_2)
  CI_power
  
  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_beta_reestimate[iter] = betahat
  iter_power_reestimate[iter] = 1 - betahat
  
  iter_CI_beta_low[iter] = CI_beta[1]
  iter_CI_beta_high[iter] = CI_beta[2]
  iter_CI_power_low[iter] = CI_power[1]
  iter_CI_power_high[iter] = CI_power[2]
  
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_betahat = 1 - mean(non_ssr_rejectnull)
  
  non_ssr_standard_error_1 = sqrt((non_ssr_betahat*(1 - non_ssr_betahat))/Nsim)
  non_ssr_standard_error_1
  
  non_ssr_CI_beta = c(non_ssr_betahat - non_ssr_standard_error_1, non_ssr_betahat + non_ssr_standard_error_1)
  
  non_ssr_iter_beta_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_beta_low[iter] = non_ssr_CI_beta[1]
  non_ssr_iter_CI_beta_high[iter] = non_ssr_CI_beta[2]
  
  
  # non ssr power
  non_ssr_power = 1 - non_ssr_betahat
  
  non_ssr_standard_error_2 = sqrt((non_ssr_power*(1 - non_ssr_power))/Nsim)
  non_ssr_standard_error_2
  
  
  # Confidence Intervals
  non_ssr_CI_power = c(non_ssr_power - non_ssr_standard_error_2, non_ssr_power + non_ssr_standard_error_2)
  non_ssr_CI_power
  
  
  non_ssr_iter_power_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_power_low[iter] = non_ssr_CI_power[1]
  non_ssr_iter_CI_power_high[iter] = non_ssr_CI_power[2]
  
}


results_h1_k0.25 = data.frame("delta" = rep(2.34,rows),
                              "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                              "Sigma Re-estimate" = iter_sigma_reestimate, 
                              "Sample Size Initial" = non_ssr_ss, "Sample Size Re-estimate" = iter_ss_reestimate, 
                              "Beta Re-estimate" = iter_beta_reestimate, 
                              "Confidence Interval Beta Low" = iter_CI_beta_low, "Confidence Interval Beta High" = iter_CI_beta_high,
                              "NON SSR Beta Re-estimate" = non_ssr_iter_beta_reestimate, 
                              "NON SSR Confidence Interval Beta Low" = non_ssr_iter_CI_beta_low, "NON SSR Confidence Interval Beta High" = non_ssr_iter_CI_beta_low,
                              "Power Re-estimate" = iter_power_reestimate,
                              "Confidence Interval Power Low" = iter_CI_power_low, "Confidence Interval Power High" = iter_CI_power_high,
                              "Power Re-estimate" = non_ssr_iter_power_reestimate,
                              "NON SSR Confidence Interval Power Low" = non_ssr_iter_CI_power_low, "NON SSR Confidence Interval Power High" = non_ssr_iter_CI_power_high)
print(results_h1_k0.25)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_h1_k0.25$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_h1_k0.25$Sample.Size.Re.estimate = ss_reestimate


write.csv(results_h1_k0.25, "Results/allocation_study_blinded_h1_k_0;25.csv", row.names = FALSE)

write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h1_k_0;25_distribution.csv", row.names = FALSE)














### SECTION II: SSR - Chow, Wang and Shao - k = 0.50 - (where H1 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H1 holds



iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )    

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_beta_reestimate = rep(NA, rows)
iter_power_reestimate = rep(NA, rows)


iter_CI_beta_low = rep(NA, rows)
iter_CI_beta_high = rep(NA, rows)
iter_CI_power_low = rep(NA, rows)
iter_CI_power_high = rep(NA, rows)


df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_beta_reestimate = rep(NA, rows)
non_ssr_iter_power_reestimate = rep(NA, rows)

non_ssr_iter_CI_beta_low  = rep(NA, rows)
non_ssr_iter_CI_beta_high = rep(NA, rows)

non_ssr_iter_CI_power_low  = rep(NA, rows)
non_ssr_iter_CI_power_high = rep(NA, rows)


for(iter in 1:rows)
{

  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters 
  sigma = iter_sigma[iter]
  
  
  Delta = 2.34
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 2.34 # The ACTUAL treatment effect, under the alternate hypothesis.
  
  

  k = 0.50  # Allocation ratio
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.  
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) 
  sigma_reestimates = rep(NA, Nsim) 
  N_reestimates = rep(NA, Nsim)
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    
    
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                     
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
  betahat = 1 - mean(rejectnull)
  betahat

  standard_error_1 = sqrt((betahat*(1 - betahat))/Nsim)
  standard_error_1
  
  
  # Confidence Intervals
  CI_beta = c(betahat - standard_error_1, betahat + standard_error_1)
  CI_beta
  
  
  sigma_reestimate = mean(sigma_reestimates)
  sigma_reestimate

  power = 1 - betahat
  cat("Power : ",power)
  
  standard_error_2 = sqrt((power*(1 - power))/Nsim)
  standard_error_2
  
  
  # Confidence Intervals
  CI_power = c(power - standard_error_2, power + standard_error_2)
  CI_power
  
  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_beta_reestimate[iter] = betahat
  iter_power_reestimate[iter] = 1 - betahat
  
  iter_CI_beta_low[iter] = CI_beta[1]
  iter_CI_beta_high[iter] = CI_beta[2]
  iter_CI_power_low[iter] = CI_power[1]
  iter_CI_power_high[iter] = CI_power[2]
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_betahat = 1 - mean(non_ssr_rejectnull)
  
  non_ssr_standard_error_1 = sqrt((non_ssr_betahat*(1 - non_ssr_betahat))/Nsim)
  non_ssr_standard_error_1
  
  non_ssr_CI_beta = c(non_ssr_betahat - non_ssr_standard_error_1, non_ssr_betahat + non_ssr_standard_error_1)
  
  non_ssr_iter_beta_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_beta_low[iter] = non_ssr_CI_beta[1]
  non_ssr_iter_CI_beta_high[iter] = non_ssr_CI_beta[2]
  
  
  # non ssr power
  non_ssr_power = 1 - non_ssr_betahat
  
  non_ssr_standard_error_2 = sqrt((non_ssr_power*(1 - non_ssr_power))/Nsim)
  non_ssr_standard_error_2
  
  
  # Confidence Intervals
  non_ssr_CI_power = c(non_ssr_power - non_ssr_standard_error_2, non_ssr_power + non_ssr_standard_error_2)
  non_ssr_CI_power
  
  
  non_ssr_iter_power_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_power_low[iter] = non_ssr_CI_power[1]
  non_ssr_iter_CI_power_high[iter] = non_ssr_CI_power[2]
  
  
}


results_h1_k0.50 = data.frame("delta" = rep(2.34,rows),
                              "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                              "Sigma Re-estimate" = iter_sigma_reestimate, 
                              "Sample Size Initial" = non_ssr_ss, "Sample Size Re-estimate" = iter_ss_reestimate, 
                              "Beta Re-estimate" = iter_beta_reestimate, 
                              "Confidence Interval Beta Low" = iter_CI_beta_low, "Confidence Interval Beta High" = iter_CI_beta_high,
                              "NON SSR Beta Re-estimate" = non_ssr_iter_beta_reestimate, 
                              "NON SSR Confidence Interval Beta Low" = non_ssr_iter_CI_beta_low, "NON SSR Confidence Interval Beta High" = non_ssr_iter_CI_beta_low,
                              "Power Re-estimate" = iter_power_reestimate,
                              "Confidence Interval Power Low" = iter_CI_power_low, "Confidence Interval Power High" = iter_CI_power_high,
                              "Power Re-estimate" = non_ssr_iter_power_reestimate,
                              "NON SSR Confidence Interval Power Low" = non_ssr_iter_CI_power_low, "NON SSR Confidence Interval Power High" = non_ssr_iter_CI_power_high)
print(results_h1_k0.50)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_h1_k0.50$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_h1_k0.50$Sample.Size.Re.estimate = ss_reestimate



write.csv(results_h1_k0.50, "Results/allocation_study_blinded_h1_k_0;50.csv", row.names = FALSE)

write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h1_k_0;50_distribution.csv", row.names = FALSE)













### SECTION III: SSR - Chow, Wang and Shao - k = 0.75 - (where H1 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H1 holds




iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )   

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_beta_reestimate = rep(NA, rows)
iter_power_reestimate = rep(NA, rows)


iter_CI_beta_low = rep(NA, rows)
iter_CI_beta_high = rep(NA, rows)
iter_CI_power_low = rep(NA, rows)
iter_CI_power_high = rep(NA, rows)

df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_beta_reestimate = rep(NA, rows)
non_ssr_iter_power_reestimate = rep(NA, rows)

non_ssr_iter_CI_beta_low  = rep(NA, rows)
non_ssr_iter_CI_beta_high = rep(NA, rows)

non_ssr_iter_CI_power_low  = rep(NA, rows)
non_ssr_iter_CI_power_high = rep(NA, rows)

for(iter in 1:rows)
{

  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters
  sigma = iter_sigma[iter]
  
  
  Delta = 2.34
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 2.34 # The ACTUAL treatment effect, under the alternate hypothesis.
  
  

  k = 0.75  # Allocation ratio
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   
  
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) 
  sigma_reestimates = rep(NA, Nsim) 
  N_reestimates = rep(NA, Nsim)
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    
    
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                     
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    # t-test on re-estimated sample size data
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
  betahat = 1 - mean(rejectnull)
  betahat

  standard_error_1 = sqrt((betahat*(1 - betahat))/Nsim)
  standard_error_1
  
  
  # Confidence Intervals
  CI_beta = c(betahat - standard_error_1, betahat + standard_error_1)
  CI_beta
  
  
  sigma_reestimate = mean(sigma_reestimates)
  sigma_reestimate

  
  power = 1 - betahat
  cat("Power : ",power)
  
  standard_error_2 = sqrt((power*(1 - power))/Nsim)
  standard_error_2
  
  
  # Confidence Intervals
  CI_power = c(power - standard_error_2, power + standard_error_2)
  CI_power
  
  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_beta_reestimate[iter] = betahat
  iter_power_reestimate[iter] = 1 - betahat
  
  iter_CI_beta_low[iter] = CI_beta[1]
  iter_CI_beta_high[iter] = CI_beta[2]
  iter_CI_power_low[iter] = CI_power[1]
  iter_CI_power_high[iter] = CI_power[2]
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_betahat = 1 - mean(non_ssr_rejectnull)
  
  non_ssr_standard_error_1 = sqrt((non_ssr_betahat*(1 - non_ssr_betahat))/Nsim)
  non_ssr_standard_error_1
  
  non_ssr_CI_beta = c(non_ssr_betahat - non_ssr_standard_error_1, non_ssr_betahat + non_ssr_standard_error_1)
  
  non_ssr_iter_beta_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_beta_low[iter] = non_ssr_CI_beta[1]
  non_ssr_iter_CI_beta_high[iter] = non_ssr_CI_beta[2]
  
  
  # non ssr power
  non_ssr_power = 1 - non_ssr_betahat
  
  non_ssr_standard_error_2 = sqrt((non_ssr_power*(1 - non_ssr_power))/Nsim)
  non_ssr_standard_error_2
  
  
  # Confidence Intervals
  non_ssr_CI_power = c(non_ssr_power - non_ssr_standard_error_2, non_ssr_power + non_ssr_standard_error_2)
  non_ssr_CI_power
  
  
  non_ssr_iter_power_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_power_low[iter] = non_ssr_CI_power[1]
  non_ssr_iter_CI_power_high[iter] = non_ssr_CI_power[2]
  
  
  
}


results_h1_k0.75 = data.frame("delta" = rep(2.34,rows),
                              "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                              "Sigma Re-estimate" = iter_sigma_reestimate, 
                              "Sample Size Initial" = non_ssr_ss, "Sample Size Re-estimate" = iter_ss_reestimate, 
                              "Beta Re-estimate" = iter_beta_reestimate, 
                              "Confidence Interval Beta Low" = iter_CI_beta_low, "Confidence Interval Beta High" = iter_CI_beta_high,
                              "NON SSR Beta Re-estimate" = non_ssr_iter_beta_reestimate, 
                              "NON SSR Confidence Interval Beta Low" = non_ssr_iter_CI_beta_low, "NON SSR Confidence Interval Beta High" = non_ssr_iter_CI_beta_low,
                              "Power Re-estimate" = iter_power_reestimate,
                              "Confidence Interval Power Low" = iter_CI_power_low, "Confidence Interval Power High" = iter_CI_power_high,
                              "Power Re-estimate" = non_ssr_iter_power_reestimate,
                              "NON SSR Confidence Interval Power Low" = non_ssr_iter_CI_power_low, "NON SSR Confidence Interval Power High" = non_ssr_iter_CI_power_high)
print(results_h1_k0.75)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_h1_k0.75$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_h1_k0.75$Sample.Size.Re.estimate = ss_reestimate



write.csv(results_h1_k0.75, "Results/allocation_study_blinded_h1_k_0;75.csv", row.names = FALSE)

write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h1_k_0;75_distribution.csv", row.names = FALSE)













### SECTION IV: SSR - Chow, Wang and Shao - k = 1.0 - (where H1 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H1 holds



 

iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )   

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_beta_reestimate = rep(NA, rows)
iter_power_reestimate = rep(NA, rows)


iter_CI_beta_low = rep(NA, rows)
iter_CI_beta_high = rep(NA, rows)
iter_CI_power_low = rep(NA, rows)
iter_CI_power_high = rep(NA, rows)


df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_beta_reestimate = rep(NA, rows)
non_ssr_iter_power_reestimate = rep(NA, rows)

non_ssr_iter_CI_beta_low  = rep(NA, rows)
non_ssr_iter_CI_beta_high = rep(NA, rows)

non_ssr_iter_CI_power_low  = rep(NA, rows)
non_ssr_iter_CI_power_high = rep(NA, rows)


for(iter in 1:rows)
{
  
  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters
  sigma = iter_sigma[iter]
  
  
  Delta = 2.34 
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 2.34 # The ACTUAL treatment effect, under the alternate hypothesis.
  
  

  k = 1  # Allocation ratio 
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   
  
  
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim)
  sigma_reestimates = rep(NA, Nsim) 
  N_reestimates = rep(NA, Nsim) 
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
 
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    
    
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                     
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    # t-test on re-estimated sample size data
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
  betahat = 1 - mean(rejectnull)
  betahat

  standard_error_1 = sqrt((betahat*(1 - betahat))/Nsim)
  standard_error_1
  
  
  # Confidence Intervals
  CI_beta = c(betahat - standard_error_1, betahat + standard_error_1)
  CI_beta
  
  
  sigma_reestimate = mean(sigma_reestimates)
  sigma_reestimate
  # 500000 simulations: 6.479028
  
  power = 1 - betahat
  cat("Power : ",power)
  
  standard_error_2 = sqrt((power*(1 - power))/Nsim)
  standard_error_2
  
  
  # Confidence Intervals
  CI_power = c(power - standard_error_2, power + standard_error_2)
  CI_power
  
  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_beta_reestimate[iter] = betahat
  iter_power_reestimate[iter] = 1 - betahat
  
  iter_CI_beta_low[iter] = CI_beta[1]
  iter_CI_beta_high[iter] = CI_beta[2]
  iter_CI_power_low[iter] = CI_power[1]
  iter_CI_power_high[iter] = CI_power[2]
  
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_betahat = 1 - mean(non_ssr_rejectnull)
  
  non_ssr_standard_error_1 = sqrt((non_ssr_betahat*(1 - non_ssr_betahat))/Nsim)
  non_ssr_standard_error_1
  
  non_ssr_CI_beta = c(non_ssr_betahat - non_ssr_standard_error_1, non_ssr_betahat + non_ssr_standard_error_1)
  
  non_ssr_iter_beta_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_beta_low[iter] = non_ssr_CI_beta[1]
  non_ssr_iter_CI_beta_high[iter] = non_ssr_CI_beta[2]
  
  
  # non ssr power
  non_ssr_power = 1 - non_ssr_betahat
  
  non_ssr_standard_error_2 = sqrt((non_ssr_power*(1 - non_ssr_power))/Nsim)
  non_ssr_standard_error_2
  
  
  # Confidence Intervals
  non_ssr_CI_power = c(non_ssr_power - non_ssr_standard_error_2, non_ssr_power + non_ssr_standard_error_2)
  non_ssr_CI_power
  
  
  non_ssr_iter_power_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_power_low[iter] = non_ssr_CI_power[1]
  non_ssr_iter_CI_power_high[iter] = non_ssr_CI_power[2]
  
  
}


results_h1_k1.00 = data.frame("delta" = rep(2.34,rows),
                              "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                              "Sigma Re-estimate" = iter_sigma_reestimate, 
                              "Sample Size Initial" = non_ssr_ss, "Sample Size Re-estimate" = iter_ss_reestimate, 
                              "Beta Re-estimate" = iter_beta_reestimate, 
                              "Confidence Interval Beta Low" = iter_CI_beta_low, "Confidence Interval Beta High" = iter_CI_beta_high,
                              "NON SSR Beta Re-estimate" = non_ssr_iter_beta_reestimate, 
                              "NON SSR Confidence Interval Beta Low" = non_ssr_iter_CI_beta_low, "NON SSR Confidence Interval Beta High" = non_ssr_iter_CI_beta_low,
                              "Power Re-estimate" = iter_power_reestimate,
                              "Confidence Interval Power Low" = iter_CI_power_low, "Confidence Interval Power High" = iter_CI_power_high,
                              "Power Re-estimate" = non_ssr_iter_power_reestimate,
                              "NON SSR Confidence Interval Power Low" = non_ssr_iter_CI_power_low, "NON SSR Confidence Interval Power High" = non_ssr_iter_CI_power_high)
print(results_h1_k1.00)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_h1_k1.00$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_h1_k1.00$Sample.Size.Re.estimate = ss_reestimate



write.csv(results_h1_k1.00, "Results/allocation_study_blinded_h1_k_1;00.csv", row.names = FALSE)

write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h1_k_1;00_distribution.csv", row.names = FALSE)










### SECTION V: SSR - Chow, Wang and Shao - k = 1.5 - (where H1 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H1 holds

# set.seed(048) # Set a seed to ensure reproducible results.




iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )    

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_beta_reestimate = rep(NA, rows)
iter_power_reestimate = rep(NA, rows)


iter_CI_beta_low = rep(NA, rows)
iter_CI_beta_high = rep(NA, rows)
iter_CI_power_low = rep(NA, rows)
iter_CI_power_high = rep(NA, rows)


df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_beta_reestimate = rep(NA, rows)
non_ssr_iter_power_reestimate = rep(NA, rows)

non_ssr_iter_CI_beta_low  = rep(NA, rows)
non_ssr_iter_CI_beta_high = rep(NA, rows)

non_ssr_iter_CI_power_low  = rep(NA, rows)
non_ssr_iter_CI_power_high = rep(NA, rows)


for(iter in 1:rows)
{
   
  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters
  sigma = iter_sigma[iter]
  
  
  Delta = 2.34
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 2.34 # The ACTUAL treatment effect, under the alternate hypothesis.
  
  

  k = 1.5  # Allocation ratio 
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.  
  
  
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) 
  sigma_reestimates = rep(NA, Nsim)
  N_reestimates = rep(NA, Nsim) 
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  # Loop to perform our simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
  
    
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    
    
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                     
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
  betahat = 1 - mean(rejectnull)
  betahat
  
  
  standard_error_1 = sqrt((betahat*(1 - betahat))/Nsim)
  standard_error_1
  
  
  # Confidence Intervals
  CI_beta = c(betahat - standard_error_1, betahat + standard_error_1)
  CI_beta
  
  
  sigma_reestimate = mean(sigma_reestimates)
  sigma_reestimate

  
  power = 1 - betahat
  cat("Power : ",power)
  
  standard_error_2 = sqrt((power*(1 - power))/Nsim)
  standard_error_2
  
  
  # Confidence Intervals
  CI_power = c(power - standard_error_2, power + standard_error_2)
  CI_power
  
  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_beta_reestimate[iter] = betahat
  iter_power_reestimate[iter] = 1 - betahat
  
  iter_CI_beta_low[iter] = CI_beta[1]
  iter_CI_beta_high[iter] = CI_beta[2]
  iter_CI_power_low[iter] = CI_power[1]
  iter_CI_power_high[iter] = CI_power[2]
  
  
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_betahat = 1 - mean(non_ssr_rejectnull)
  
  non_ssr_standard_error_1 = sqrt((non_ssr_betahat*(1 - non_ssr_betahat))/Nsim)
  non_ssr_standard_error_1
  
  non_ssr_CI_beta = c(non_ssr_betahat - non_ssr_standard_error_1, non_ssr_betahat + non_ssr_standard_error_1)
  
  non_ssr_iter_beta_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_beta_low[iter] = non_ssr_CI_beta[1]
  non_ssr_iter_CI_beta_high[iter] = non_ssr_CI_beta[2]
  
  
  # non ssr power
  non_ssr_power = 1 - non_ssr_betahat
  
  non_ssr_standard_error_2 = sqrt((non_ssr_power*(1 - non_ssr_power))/Nsim)
  non_ssr_standard_error_2
  
  
  # Confidence Intervals
  non_ssr_CI_power = c(non_ssr_power - non_ssr_standard_error_2, non_ssr_power + non_ssr_standard_error_2)
  non_ssr_CI_power
  
  
  non_ssr_iter_power_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_power_low[iter] = non_ssr_CI_power[1]
  non_ssr_iter_CI_power_high[iter] = non_ssr_CI_power[2]
  
  
}


results_h1_k1.50 = data.frame("delta" = rep(2.34,rows),
                              "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                              "Sigma Re-estimate" = iter_sigma_reestimate, 
                              "Sample Size Initial" = non_ssr_ss, "Sample Size Re-estimate" = iter_ss_reestimate, 
                              "Beta Re-estimate" = iter_beta_reestimate, 
                              "Confidence Interval Beta Low" = iter_CI_beta_low, "Confidence Interval Beta High" = iter_CI_beta_high,
                              "NON SSR Beta Re-estimate" = non_ssr_iter_beta_reestimate, 
                              "NON SSR Confidence Interval Beta Low" = non_ssr_iter_CI_beta_low, "NON SSR Confidence Interval Beta High" = non_ssr_iter_CI_beta_low,
                              "Power Re-estimate" = iter_power_reestimate,
                              "Confidence Interval Power Low" = iter_CI_power_low, "Confidence Interval Power High" = iter_CI_power_high,
                              "Power Re-estimate" = non_ssr_iter_power_reestimate,
                              "NON SSR Confidence Interval Power Low" = non_ssr_iter_CI_power_low, "NON SSR Confidence Interval Power High" = non_ssr_iter_CI_power_high)
print(results_h1_k1.50)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_h1_k1.50$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_h1_k1.50$Sample.Size.Re.estimate = ss_reestimate



write.csv(results_h1_k1.50, "Results/allocation_study_blinded_h1_k_1;50.csv", row.names = FALSE)

write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h1_k_1;50_distribution.csv", row.names = FALSE)





















### SECTION VI: SSR - Chow, Wang and Shao - k = 2.0 - (where H1 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H1 holds



 

iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0 )    

rows = length(iter_sigma)
iter_sigma_reestimate = rep(NA,rows)
iter_ss_reestimate = rep(NA, rows)
iter_beta_reestimate = rep(NA, rows)
iter_power_reestimate = rep(NA, rows)


iter_CI_beta_low = rep(NA, rows)
iter_CI_beta_high = rep(NA, rows)
iter_CI_power_low = rep(NA, rows)
iter_CI_power_high = rep(NA, rows)


df_iter_ss_reestimate = data.frame("Iteration" = 1:Nsim)
non_ssr_ss = rep(NA,rows)
non_ssr_iter_beta_reestimate = rep(NA, rows)
non_ssr_iter_power_reestimate = rep(NA, rows)

non_ssr_iter_CI_beta_low  = rep(NA, rows)
non_ssr_iter_CI_beta_high = rep(NA, rows)

non_ssr_iter_CI_power_low  = rep(NA, rows)
non_ssr_iter_CI_power_high = rep(NA, rows)


for(iter in 1:rows)
{

  sigmaTRUE = iter_sigmaTRUE[iter]
  
  
  # Assumed parameters
  sigma = iter_sigma[iter]
  
  
  Delta = 2.34 
  
  alpha = 0.05 # Desired error rates for our
  beta = 0.2   # simulated trial.
  
  delta = 2.34 # The ACTUAL treatment effect, under the alternate hypothesis.
  
  

  k = 2  # Allocation ratio
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = round(k * n2)  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.
  
  
  # Number of simulated trials to run.
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) 
  sigma_reestimates = rep(NA, Nsim) 
  N_reestimates = rep(NA, Nsim)
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  # Loop to perform our simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim1 = round(n1/2, 0)                
    interim2 = round(n2/2, 0)
    
  
    interim_samples = c(resultsC[1:interim1], resultsT[1:interim2])
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = round(k * n2_reestimate) # n1 is control size
    
    
    N_reestimate = n1_reestimate + n2_reestimate
    N_reestimates[i] = N_reestimate
    cat("\n\nTrial:",i,"\t\tRe-estimated Sample Size: ",N_reestimate)
    
    
    
    
    if(sigma_reestimate > sigma)   # We under-estimated sigma at the start (SSR will increase ss)
    {
      if(n1_reestimate > n1)
      { 
        difference1 = n1_reestimate - n1
      }
      else
      { 
        difference1 = 0
      }                                    
      
      if(n2_reestimate > n2)
      { 
        difference2 = n2_reestimate - n2 
      }
      else
      { 
        difference2 = 0 
      } 
      
      resultsC_extension = rnorm(difference1, muC, sigmaTRUE)
      resultsT_extension = rnorm(difference2, muT, sigmaTRUE)
      
      resultsC_total = c(resultsC, resultsC_extension)
      resultsT_total = c(resultsT, resultsT_extension)
      
    }
    
    if(sigma_reestimate < sigma)   # We over-estimated sigma at the start (SSR will decrease ss)
    {
      
      if(n1_reestimate < n1)
      { 
        difference1 = n1 - n1_reestimate 
        resultsC_total = head(resultsC, -difference1) 
      }
      else
      { 
        difference1 = 0
        resultsC_total = resultsC 
      }                                     
      
      if(n2_reestimate < n2)
      { 
        difference2 = n2 - n2_reestimate 
        resultsT_total = head(resultsT, -difference2) 
      }
      else
      { 
        difference2 = 0 
        resultsT_total = resultsT 
      } 
      
      
      
    }
    
    outcome = t.test(resultsC_total, resultsT_total, var.equal = TRUE, alternative = "greater")
    pvalues[i] = outcome$p.value
    
    outcome$p.value < 0.05
    
    
    # t-test on original data (without re-estimation)
    non_ssr_outcome = t.test(resultsC, resultsT, var.equal = TRUE, alternative = "greater")
    non_ssr_pvalues[i] = non_ssr_outcome$p.value
    
    
    cat("\nIteration end")
    
  }
  
  
  
  rejectnull = (pvalues < 0.05) # Vector of conclusions in each simulated trial.
  
  # To determine the proportion of these trials that INCORRECTLY
  # reject the null hypothesis, we take the arithmetic mean:
  
  betahat = 1 - mean(rejectnull)
  betahat

  standard_error_1 = sqrt((betahat*(1 - betahat))/Nsim)
  standard_error_1
  
  
  # Confidence Intervals
  CI_beta = c(betahat - standard_error_1, betahat + standard_error_1)
  CI_beta
  
  
  sigma_reestimate = mean(sigma_reestimates)
  sigma_reestimate

  power = 1 - betahat
  cat("Power : ",power)
  
  standard_error_2 = sqrt((power*(1 - power))/Nsim)
  standard_error_2
  
  
  # Confidence Intervals
  CI_power = c(power - standard_error_2, power + standard_error_2)
  CI_power
  
  
  cat("Re-estimated Total Sample Size :", round(mean(N_reestimates)))
  
  
  
  
  
  iter_sigma_reestimate[iter] = sigma_reestimate
  iter_ss_reestimate[iter] = round(mean(N_reestimates))
  iter_beta_reestimate[iter] = betahat
  iter_power_reestimate[iter] = 1 - betahat
  
  iter_CI_beta_low[iter] = CI_beta[1]
  iter_CI_beta_high[iter] = CI_beta[2]
  iter_CI_power_low[iter] = CI_power[1]
  iter_CI_power_high[iter] = CI_power[2]
  
  
 
  # Collecting the re-estimated sample size to visualize the distribution 
  df_iter_ss_reestimate[paste("Combination_",iter, sep="")] = N_reestimates
  
  non_ssr_ss[iter] = N
  
  non_ssr_rejectnull = (non_ssr_pvalues < 0.05)
  non_ssr_betahat = 1 - mean(non_ssr_rejectnull)
  
  non_ssr_standard_error_1 = sqrt((non_ssr_betahat*(1 - non_ssr_betahat))/Nsim)
  non_ssr_standard_error_1
  
  non_ssr_CI_beta = c(non_ssr_betahat - non_ssr_standard_error_1, non_ssr_betahat + non_ssr_standard_error_1)
  
  non_ssr_iter_beta_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_beta_low[iter] = non_ssr_CI_beta[1]
  non_ssr_iter_CI_beta_high[iter] = non_ssr_CI_beta[2]
  
  
  # non ssr power
  non_ssr_power = 1 - non_ssr_betahat
  
  non_ssr_standard_error_2 = sqrt((non_ssr_power*(1 - non_ssr_power))/Nsim)
  non_ssr_standard_error_2
  
  
  # Confidence Intervals
  non_ssr_CI_power = c(non_ssr_power - non_ssr_standard_error_2, non_ssr_power + non_ssr_standard_error_2)
  non_ssr_CI_power
  
  
  non_ssr_iter_power_reestimate[iter] = non_ssr_betahat
  
  non_ssr_iter_CI_power_low[iter] = non_ssr_CI_power[1]
  non_ssr_iter_CI_power_high[iter] = non_ssr_CI_power[2]
  
}


results_h1_k2.00 = data.frame("delta" = rep(2.34,rows),
                              "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                              "Sigma Re-estimate" = iter_sigma_reestimate, 
                              "Sample Size Initial" = non_ssr_ss, "Sample Size Re-estimate" = iter_ss_reestimate, 
                              "Beta Re-estimate" = iter_beta_reestimate, 
                              "Confidence Interval Beta Low" = iter_CI_beta_low, "Confidence Interval Beta High" = iter_CI_beta_high,
                              "NON SSR Beta Re-estimate" = non_ssr_iter_beta_reestimate, 
                              "NON SSR Confidence Interval Beta Low" = non_ssr_iter_CI_beta_low, "NON SSR Confidence Interval Beta High" = non_ssr_iter_CI_beta_low,
                              "Power Re-estimate" = iter_power_reestimate,
                              "Confidence Interval Power Low" = iter_CI_power_low, "Confidence Interval Power High" = iter_CI_power_high,
                              "Power Re-estimate" = non_ssr_iter_power_reestimate,
                              "NON SSR Confidence Interval Power Low" = non_ssr_iter_CI_power_low, "NON SSR Confidence Interval Power High" = non_ssr_iter_CI_power_high)
print(results_h1_k2.00)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,rows)
for(i in 1:rows){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results_h1_k2.00$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results_h1_k2.00$Sample.Size.Re.estimate = ss_reestimate



write.csv(results_h1_k2.00, "Results/allocation_study_blinded_h1_k_2;00.csv", row.names = FALSE)

write.csv(df_iter_ss_reestimate, "Results/allocation_study_blinded_h1_k_2;00_distribution.csv", row.names = FALSE)









