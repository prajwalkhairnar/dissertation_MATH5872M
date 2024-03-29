# METHOD 1

# Chow S., Lokhnygina Y., Shao J. and Wang H. (2017). Sample size calculations in clinical research 
# (Third edition.). Taylor and Francis


# set.seed(048) # Set a seed to ensure reproducible results.



### SECTION I: SSR - Chow, Wang and Shao - k = 1 - (where H0 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H0 holds


iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5, 5.0, 5.0, 5.0, 5.0, 5.0 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0, 4.5, 4.8, 5.0, 5.3, 5.5 )   

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
  beta = 0.2   # SIMULATED trial.
  
  delta = 0 # The ACTUAL treatment effect, under the null hypothesis.
  
  
  # Calculate the desired sample size using the formula (for one group): 
  
  k = 1  # Allocation ratio 
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = k * n2  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment. 
  
  
  Nsim = 500000
  
  
  
  pvalues = rep(NA, Nsim) # An empty vector to store the p-values.
  sigma_reestimates = rep(NA, Nsim) # An empty vector to store the re-estimated sigma values.
  N_reestimates = rep(NA, Nsim) # N re-estimated values
  
  
  non_ssr_pvalues = rep(NA, Nsim)
  
  
  # Loop to perform our simulations:
  for(i in 1:Nsim){
    
    resultsC = rnorm(n1, muC, sigmaTRUE)
    resultsT = rnorm(n2, muT, sigmaTRUE)
    interim = round(n2/2, 0)                # n1 and n2 are same
    
    
    
    
    interim_samples = c(resultsC[1:interim], resultsT[1:interim])  
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = k * n2_reestimate # n1 is control size
    
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
  # 500000 simulations: 6.479028
  
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



results1_1 = data.frame("delta" = rep(0,10),
                        "Sigma TRUE" = iter_sigmaTRUE, "Sigma Assumed" = iter_sigma,
                        "Sigma Re-estimate" = iter_sigma_reestimate, 
                        "Sample Size Initial" = non_ssr_ss , "Sample Size Re-estimate" = iter_ss_reestimate, 
                        "Alpha Re-estimate" = iter_alpha_reestimate,
                        "Confidence Interval Low" = iter_CI_low, "Confidence Interval High" = iter_CI_high, 
                        "NON SSR Alpha Re-estimate" = non_ssr_iter_alpha_reestimate, 
                        "NON SSR Confidence Interval Low" = non_ssr_iter_CI_low, "NON SSR Confidence Interval High" = non_ssr_iter_CI_high)

print(results1_1)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,10)
for(i in 1:10){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results1_1$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results1_1$Sample.Size.Re.estimate = ss_reestimate

write.csv(results1_1, "Results/met1_blinded_h0.csv", row.names = FALSE)


write.csv(df_iter_ss_reestimate, "Results/met1_blinded_h0_distribution.csv", row.names = FALSE)
























### SECTION II: SSR - Chow, Wang and Shao - k = 1 - (where H1 holds) - BLINDED - multiple trial runs ----


# Simulation of a trial under the assumption that H1 holds

# set.seed(048) # Set a seed to ensure reproducible results.





iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5, 5.0, 5.0, 5.0, 5.0, 5.0 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0, 4.5, 4.8, 5.0, 5.3, 5.5 )   

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
  
  
  # Calculate the desired sample size using the formula (for one group): 
  
  k = 1  # Allocation ratio (here 1:1)
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = k * n2  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.  
  
  
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
    interim = round(n2/2, 0)                # n1 and n2 are same
    
    
    
    
    
    # BLINDED Approach
    interim_samples = c(resultsC[1:interim], resultsT[1:interim])  
    
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = k * n2_reestimate # n1 is control size
    
    
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


results1_2 = data.frame("delta" = rep(2.34,10),
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
print(results1_2)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))



ss_reestimate = rep(NA,10)
for(i in 1:10){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results1_2$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results1_2$Sample.Size.Re.estimate = ss_reestimate



write.csv(results1_2, "Results/met1_blinded_h1.csv", row.names = FALSE)

write.csv(df_iter_ss_reestimate, "Results/met1_blinded_h1_distribution.csv", row.names = FALSE)


























### SECTION III: SSR - Chow, Wang and Shao - k = 1 - (where H1 holds) - Partially UNBLINDED(C) - multiple trial runs ----

# Using only control data

# Simulation of a trial under the assumption that H1 holds


iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5, 5.0, 5.0, 5.0, 5.0, 5.0 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0, 4.5, 4.8, 5.0, 5.3, 5.5 )   

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
  
  
  k = 1  # Allocation ratio 
  
  n2 = round((((qnorm(alpha) + qnorm(beta)) ^2)*(sigma^2)*(1+(1/k)))/(Delta^2))   # n2 is treatment group size
  
  n1 = k * n2  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment. 
  
  
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
    interim = round(n2/2, 0)                # n1 and n2 are same
    
    
    
    
    
    # Partially UNBLINDED - Control
    interim_samples = resultsC[1:interim]
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = k * n2_reestimate # n1 is control size
    
    
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


results1_3 = data.frame("delta" = rep(2.34,10),
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
print(results1_3)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))


ss_reestimate = rep(NA,10)
for(i in 1:10){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results1_3$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results1_3$Sample.Size.Re.estimate = ss_reestimate


write.csv(results1_3, "Results/met1_part_unblinded_C_h1.csv", row.names = FALSE)

write.csv(df_iter_ss_reestimate, "Results/met1_part_unblinded_C_h1_distribution.csv", row.names = FALSE)































### SECTION IV: SSR - Chow, Wang and Shao - k = 1 - (where H1 holds) - Partially UNBLINDED(T) - multiple trial runs ----

# Using only treatment data

# Simulation of a trial under the assumption that H1 holds





iter_sigmaTRUE = c( 6.5, 6.5, 6.5, 6.5, 6.5, 5.0, 5.0, 5.0, 5.0, 5.0 )   
iter_sigma =     c( 6.0, 6.3, 6.5, 6.8, 7.0, 4.5, 4.8, 5.0, 5.3, 5.5 )   

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
  
  n1 = k * n2  # n1 is control group size
  
  N = n1 + n2 # Total Estimated Sample Size 
  cat("\nTotal Estimated Sample Size :", N, "\n\n")
  
  muC = 11.7 # Mean hospital stay for the control (ASSUMED).
  muT = muC - delta # Mean hospital stay for the treatment.   

  
  
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
    interim = round(n2/2, 0)                # n1 and n2 are same
    
    
    
  
    # Partially UNBLINDED - Treatment 
    interim_samples = resultsT[1:interim]
    
    
    sigma_reestimate = sqrt(var(interim_samples))
    sigma_reestimates[i] = sigma_reestimate
    
    
    n2_reestimate = round((((qnorm(alpha) + qnorm(beta))^2)*(sigma_reestimate^2)*(1+(1/k)))/(Delta^2))  # n2 is treatment group size
    n1_reestimate = k * n2_reestimate # n1 is control size
    
    
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



results1_4 = data.frame("delta" = rep(2.34,10),
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
print(results1_4)
df_iter_ss_reestimate = subset(df_iter_ss_reestimate, select = -c(Iteration))





ss_reestimate = rep(NA,10)
for(i in 1:10){
  ss_std_error = sqrt(var(df_iter_ss_reestimate[i])/Nsim)
  print(ss_std_error)
  ss_CI = c(round(mean(df_iter_ss_reestimate[,i]))- ss_std_error, 
            round(mean(df_iter_ss_reestimate[,i]))+ ss_std_error)
  print(ss_CI)
  ss_reestimate[i] = paste(results1_4$Sample.Size.Re.estimate[i], " (" , round(ss_CI[1],3), ",",  round(ss_CI[2],3), ")", sep = "")
}

results1_4$Sample.Size.Re.estimate = ss_reestimate


write.csv(results1_4, "Results/met1_part_unblinded_T_h1.csv", row.names = FALSE)



write.csv(df_iter_ss_reestimate, "Results/met1_part_unblinded_T_h1_distribution.csv", row.names = FALSE)



