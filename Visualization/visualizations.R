library(ggplot2)
library(xtable)


## Confidence Interval plot function ----


plot_conf_intrv = function(estimate, upper, lower, footnote, flag)
{
  analysis = c("sigmaTRUE = 6.5 | sigmaAssumed = 6.0",
               "sigmaTRUE = 6.5 | sigmaAssumed = 6.3",
               "sigmaTRUE = 6.5 | sigmaAssumed = 6.5",
               "sigmaTRUE = 6.5 | sigmaAssumed = 6.8",
               "sigmaTRUE = 6.5 | sigmaAssumed = 7.0",
               "sigmaTRUE = 5.0 | sigmaAssumed = 4.5",
               "sigmaTRUE = 5.0 | sigmaAssumed = 4.8",
               "sigmaTRUE = 5.0 | sigmaAssumed = 5.0",
               "sigmaTRUE = 5.0 | sigmaAssumed = 5.3",
               "sigmaTRUE = 5.0 | sigmaAssumed = 5.5"
  )
  
  
  # analysis= c("M = 10",        # Used only for Method 2 : Analysis of parameter M
  #             "M = 20",
  #             "M = 30",
  #             "M = 40",
  #             "M = 50",
  #             "M = 100",
  #             "M = 120",
  #             "M = 150")
  # 
  
  # analysis = c("k = 0.25",       # Used only for allocation ratio
  #              "k = 0.50",
  #              "k = 0.75",
  #              "k = 1.00",
  #              "k = 1.50",
  #              "k = 2.00")

  
  
  par(mar = c(6,10,1,10))    
  
  if(flag == 0){ 
    xlims = c(0.045, 0.055)
    guides = c(0.048, 0.05, 0.052)
  }
  else if(flag == 1){ 
    xlims = c(0.18, 0.22)
    guides = c(0.19, 0.20, 0.21)
  }
  else if(flag == 2){ 
    xlims = c(0.795, 0.805)
    guides = c(0.798, 0.80, 0.802)
  }
  
  
  plot(x = 0,                                  # One point at (0,0).
       xlim = xlims, ylim=c(0, 10),            # Axis limits.
       type = "n", xaxt = "n", yaxt="n",       # No points, no axes drawn.
       xlab = NULL, ylab= NULL, ann = FALSE,   # No axis labels or numbers.
       bty="n")                                # No box.
  axis(side = 1, cex.axis = 1) 
  mtext(footnote, 
        side = 1, line = 4)
  
  
  
  for(i in guides){
    
    lines(c(i, i), c(0, 11), lty = 2, col = "gray53")
    
  } 
  
  verticalpos = rev(1:10)
  # verticalpos = rev(1:8)                      # Used only for Method 2 : Analysis of parameter M
  # verticalpos = rev(1:6)                        # Used only for Allocation Ratio Analysis 
  
  
    
  mtext(text = analysis,  at = verticalpos, 
        side = 2, line = 9, outer = FALSE, las = 1, adj = 0)
  
  points(estimate, verticalpos, pch = 16) 
  
  for(i in rev(1:10) ){
  # for(i in rev(1:8) ){                       # Used only for Method 2 : Analysis of parameter M
  # for(i in rev(1:6) ){                       # Used only for Allocation Ratio Analysis 
    
    lines(c(lower[i], upper[i]), c(verticalpos[i], verticalpos[i]))
    
    lines(c(lower[i], lower[i]), c(verticalpos[i] + 0.2, verticalpos[i] - 0.2))
    
    lines(c(upper[i], upper[i]), c(verticalpos[i] + 0.2, verticalpos[i] - 0.2))
    
    print(i)
    
  }
  
  est <- formatC(estimate, format='f', digits = 5)
  
  
  
  L <- formatC(lower, format = 'f', digits = 5)
  U <- formatC(upper, format = 'f', digits = 5)
  
  interval <- paste("(", L, ", ", U, ")", sep = "")   # Type interval to check.
  
  
  results <- paste(est, interval)
  
  mtext(text = results, at = verticalpos, 
        side = 4, line = 9, outer = FALSE, las = 1, adj = 1)
  
  box("inner")
  
  
  
  
}





## Method 1: Base method ----


### H0 holds | blinded ----

met1_blinded_h0_v2 = read.csv("Results/met1_blinded_h0.csv")
met1_blinded_h0_distribution_v2 = read.csv("Results/met1_blinded_h0_distribution.csv")


# Scatter Plot
qplot(seq_along(met1_blinded_h0_distribution_v2$Combination_1), met1_blinded_h0_distribution_v2$Combination_1)

# Histogram
ggplot(data=met1_blinded_h0_distribution_v2, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met1_blinded_h0_distribution_v2$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 1: H0 holds", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")



# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)
footnote = "Method 1: Blinded | H0 Holds | Alpha Estimate"
plot_conf_intrv(met1_blinded_h0_v2$Alpha.Re.estimate, 
                met1_blinded_h0_v2$Confidence.Interval.High, 
                met1_blinded_h0_v2$Confidence.Interval.Low,
                footnote,
                0)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met1_blinded_h0_v2[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,0,2,2,3,0,0,5,5)), include.rownames = FALSE)





### H1 holds | blinded ----

met1_blinded_h1_v2 = read.csv("Results/met1_blinded_h1.csv")
met1_blinded_h1_distribution_v2 = read.csv("Results/met1_blinded_h1_distribution.csv")


# Scatter Plot
qplot(seq_along(met1_blinded_h1_distribution_v2$Combination_1), met1_blinded_h1_distribution_v2$Combination_1)

# Histogram
ggplot(data=met1_blinded_h1_distribution_v2, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met1_blinded_h1_distribution_v2$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 1: H1 holds", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")


# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)#
# BETA
footnote = "Method 1: Blinded | H1 Holds | Beta Estimate"
plot_conf_intrv(met1_blinded_h1_v2$Beta.Re.estimate, 
                met1_blinded_h1_v2$Confidence.Interval.Beta.High, 
                met1_blinded_h1_v2$Confidence.Interval.Beta.Low,
                footnote,
                1)
 

# POWER
footnote = "Method 1: Blinded | H1 Holds | Power Estimate"
plot_conf_intrv(met1_blinded_h1_v2$Power.Re.estimate, 
                met1_blinded_h1_v2$Confidence.Interval.Power.High, 
                met1_blinded_h1_v2$Confidence.Interval.Power.Low,
                footnote,
                2)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met1_blinded_h1_v2[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,2,2,2,3,0,0,5,5)), include.rownames = FALSE)




### H1 holds | partially unblinded (control) ----

met1_part_unblinded_C_h1_v2 = read.csv("Results/met1_part_unblinded_C_h1.csv")
met1_part_unblinded_C_h1_distribution_v2 = read.csv("Results/met1_part_unblinded_C_h1_distribution.csv")


# Scatter Plot
qplot(seq_along(met1_part_unblinded_C_h1_distribution_v2$Combination_1), met1_part_unblinded_C_h1_distribution_v2$Combination_1)

# Histogram
ggplot(data=met1_part_unblinded_C_h1_distribution_v2, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met1_part_unblinded_C_h1_distribution_v2$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 1: H1 holds (Partially Unblinded - Control)", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")


# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)#
# BETA
footnote = "Method 1: Partially Unlinded - Control | H1 Holds | Beta Estimate"
plot_conf_intrv(met1_part_unblinded_C_h1_v2$Beta.Re.estimate, 
                met1_part_unblinded_C_h1_v2$Confidence.Interval.Beta.High, 
                met1_part_unblinded_C_h1_v2$Confidence.Interval.Beta.Low,
                footnote,
                1)


# POWER
footnote = "Method 1: Blinded | H1 Holds | Power Estimate"
plot_conf_intrv(met1_part_unblinded_C_h1_v2$Power.Re.estimate, 
                met1_part_unblinded_C_h1_v2$Confidence.Interval.Power.High, 
                met1_part_unblinded_C_h1_v2$Confidence.Interval.Power.Low,
                footnote,
                2)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met1_part_unblinded_C_h1_v2[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,2,2,2,3,0,0,5,5)), include.rownames = FALSE)







### H1 holds | partially unblinded (treatment) ----

met1_part_unblinded_T_h1_v2 = read.csv("Results/met1_part_unblinded_T_h1.csv")
met1_part_unblinded_T_h1_distribution_v2 = read.csv("Results/met1_part_unblinded_T_h1_distribution.csv")


# Scatter Plot
qplot(seq_along(met1_part_unblinded_T_h1_distribution_v2$Combination_1), met1_part_unblinded_T_h1_distribution_v2$Combination_1)


# Histogram
ggplot(data=met1_part_unblinded_T_h1_distribution_v2, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met1_part_unblinded_T_h1_distribution_v2$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 1: H1 holds (Partially Unblinded - Treatment)", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")
 

# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)#
# BETA
footnote = "Method 1: Partially Unlinded - Treatment | H1 Holds | Beta Estimate"
plot_conf_intrv(met1_part_unblinded_T_h1_v2$Beta.Re.estimate, 
                met1_part_unblinded_T_h1_v2$Confidence.Interval.Beta.High, 
                met1_part_unblinded_T_h1_v2$Confidence.Interval.Beta.Low,
                footnote,
                1)


# POWER
footnote = "Method 1: Blinded | H1 Holds | Power Estimate"
plot_conf_intrv(met1_part_unblinded_T_h1_v2$Power.Re.estimate, 
                met1_part_unblinded_T_h1_v2$Confidence.Interval.Power.High, 
                met1_part_unblinded_T_h1_v2$Confidence.Interval.Power.Low,
                footnote,
                2)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met1_part_unblinded_T_h1_v2[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,2,2,2,3,0,0,5,5)), include.rownames = FALSE)











## Method 2: Betensky method ----



### Null version - studying parameter M ----


#### Type I error ----

met2_null_h0_v3 = read.csv("Results/met2_h0_null_M.csv")

met2_null_h0_v3



footnote = "Method 2: Blinded | H0 Holds | Alpha Estimate | sigmaTRUE = 6.5 | sigmaAssumed = 6.0"
plot_conf_intrv(met2_null_h0_v3$Alpha, 
                met2_null_h0_v3$Confidence.Interval.High, 
                met2_null_h0_v3$Confidence.Interval.Low,
                footnote,
                0)


#### Type II error ----



met2_null_h1_v3 = read.csv("Results/met2_h1_null_M.csv")

met2_null_h1_v3



footnote = "Method 2: Blinded | H1 Holds | Beta Estimate | sigmaTRUE = 6.5 | sigmaAssumed = 6.0"
plot_conf_intrv(met2_null_h1_v3$Beta, 
                met2_null_h1_v3$Confidence.Interval.Beta.High, 
                met2_null_h1_v3$Confidence.Interval.Beta.Low,
                footnote,
                1)




### Interim changes - v3 - Blinded H0 ----


met2_blinded_h0_v3 = read.csv("Results/met2/met2_blinded_h0_v3.csv")
met2_blinded_h0_distribution_v3 = read.csv("Results/met2/met2_blinded_h0_distribution_v3.csv")


# Scatter Plot
qplot(seq_along(met2_blinded_h0_distribution_v3$Combination_1), met2_blinded_h0_distribution_v3$Combination_1)

# Histogram
ggplot(data=met2_blinded_h0_distribution_v3, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met2_blinded_h0_distribution_v3$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 2: H0 holds", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")



# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)
footnote = "Method 2: Blinded | H0 Holds | Alpha Estimate"
plot_conf_intrv(met2_blinded_h0_v3$Alpha.Re.estimate, 
                met2_blinded_h0_v3$Confidence.Interval.High, 
                met2_blinded_h0_v3$Confidence.Interval.Low,
                footnote,
                0)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met2_blinded_h0_v3[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,0,2,2,3,0,0,5,5)), include.rownames = FALSE)




### Interim changes - v3 - Blinded H1 ----



met2_blinded_h1_v3 = read.csv("Results/met2/met2_blinded_h1_v3.csv")
met2_blinded_h1_distribution_v3 = read.csv("Results/met2/met2_blinded_h1_distribution_v3.csv")


# Scatter Plot
qplot(seq_along(met2_blinded_h1_distribution_v3$Combination_1), met2_blinded_h1_distribution_v3$Combination_1)


# Histogram
ggplot(data=met2_blinded_h1_distribution_v3, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met2_blinded_h1_distribution_v3$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 2: H1 holds (Blinded)", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")


# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)#
# BETA
footnote = "Method 2: Blinded | H1 Holds | Beta Estimate"
plot_conf_intrv(met2_blinded_h1_v3$Beta.Re.estimate, 
                met2_blinded_h1_v3$Confidence.Interval.Beta.High, 
                met2_blinded_h1_v3$Confidence.Interval.Beta.Low,
                footnote,
                1)


# POWER
footnote = "Method 4: Blinded | H1 Holds | Power Estimate"
plot_conf_intrv(met4_blinded_h1_v2$Power.Re.estimate, 
                met4_blinded_h1_v2$Confidence.Interval.Power.High, 
                met4_blinded_h1_v2$Confidence.Interval.Power.Low,
                footnote,
                2)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met2_blinded_h1_v3[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,2,2,2,3,0,0,5,5)), include.rownames = FALSE)












## Method 3: Bristol method ----



### H0 holds | blinded ----

met3_blinded_h0_v3 = read.csv("Results/met3_blinded_h0.csv")
met3_blinded_h0_distribution_v3 = read.csv("Results/met3_blinded_h0_distribution.csv")


# Scatter Plot
qplot(seq_along(met3_blinded_h0_distribution_v3$Combination_1), met3_blinded_h0_distribution_v3$Combination_1)

# Histogram
ggplot(data=met3_blinded_h0_distribution_v3, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met3_blinded_h0_distribution_v3$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 3: H0 holds", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")



# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)
footnote = "Method 3: Blinded | H0 Holds | Alpha Estimate"
plot_conf_intrv(met3_blinded_h0_v3$Alpha.Re.estimate, 
                met3_blinded_h0_v3$Confidence.Interval.High, 
                met3_blinded_h0_v3$Confidence.Interval.Low,
                footnote,
                0)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met3_blinded_h0_v3[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,0,2,2,3,0,0,5,5)), include.rownames = FALSE)


### H1 holds | blinded ----

met3_blinded_h1_v3 = read.csv("Results/met3_blinded_h1.csv")
met3_blinded_h1_distribution_v3 = read.csv("Results/met3_blinded_h1_distribution.csv")


# Scatter Plot
qplot(seq_along(met3_blinded_h1_distribution_v3$Combination_1), met3_blinded_h1_distribution_v3$Combination_1)

# Histogram
ggplot(data=met3_blinded_h1_distribution_v3, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met3_blinded_h1_distribution_v3$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 3: H1 holds (Blinded)", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")


# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)#
# BETA
footnote = "Method 3: Blinded | H1 Holds | Beta Estimate"
plot_conf_intrv(met3_blinded_h1_v3$Beta.Re.estimate, 
                met3_blinded_h1_v3$Confidence.Interval.Beta.High, 
                met3_blinded_h1_v3$Confidence.Interval.Beta.Low,
                footnote,
                1)


# POWER
footnote = "Method 1: Blinded | H1 Holds | Power Estimate"
plot_conf_intrv(met3_blinded_h1_v3$Power.Re.estimate, 
                met3_blinded_h1_v3$Confidence.Interval.Power.High, 
                met3_blinded_h1_v3$Confidence.Interval.Power.Low,
                footnote,
                2)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met3_blinded_h1_v3[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,2,2,2,3,0,0,5,5)), include.rownames = FALSE)











## Method 4: Kieser method ----

### H0 holds | blinded ----

met4_blinded_h0_v2 = read.csv("Results/met4_blinded_h0.csv")
met4_blinded_h0_distribution_v2 = read.csv("Results/met4_blinded_h0_distribution.csv")


# Scatter Plot
qplot(seq_along(met4_blinded_h0_distribution_v2$Combination_1), met4_blinded_h0_distribution_v2$Combination_1)


# Histogram
ggplot(data=met4_blinded_h0_distribution_v2, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met4_blinded_h0_distribution_v2$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 4: H0 holds", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")



# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)
footnote = "Method 4: Blinded | H0 Holds | Alpha Estimate"
plot_conf_intrv(met4_blinded_h0_v2$Alpha.Re.estimate, 
                met4_blinded_h0_v2$Confidence.Interval.High, 
                met4_blinded_h0_v2$Confidence.Interval.Low,
                footnote,
                0)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met4_blinded_h0_v2[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,0,2,2,3,0,0,5,5)), include.rownames = FALSE)




### H1 holds | blinded ----

met4_blinded_h1_v2 = read.csv("Results/met4_blinded_h1.csv")
met4_blinded_h1_distribution_v2 = read.csv("Results/met4_blinded_h1_distribution.csv")


# Scatter Plot
qplot(seq_along(met4_blinded_h1_distribution_v2$Combination_1), met4_blinded_h1_distribution_v2$Combination_1)


# Histogram
ggplot(data=met4_blinded_h1_distribution_v2, aes(Combination_1)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met4_blinded_h1_distribution_v2$Combination_1), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 4: H1 holds (Blinded)", subtitle = "sigmaTRUE = 6.5 | sigma = 6.0")+
  xlab("Re-estimated Sample Size")


# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)#
# BETA
footnote = "Method 4: Blinded | H1 Holds | Beta Estimate"
plot_conf_intrv(met4_blinded_h1_v2$Beta.Re.estimate, 
                met4_blinded_h1_v2$Confidence.Interval.Beta.High, 
                met4_blinded_h1_v2$Confidence.Interval.Beta.Low,
                footnote,
                1)


# POWER
footnote = "Method 4: Blinded | H1 Holds | Power Estimate"
plot_conf_intrv(met4_blinded_h1_v2$Power.Re.estimate, 
                met4_blinded_h1_v2$Confidence.Interval.Power.High, 
                met4_blinded_h1_v2$Confidence.Interval.Power.Low,
                footnote,
                2)



# latex table
cols = c(1,2,3,4,5,6,7,10)
print(xtable(met4_blinded_h1_v2[,cols], align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"), digits=c(0,2,2,2,3,0,0,5,5)), include.rownames = FALSE)



## ALLOCATION RATIO ----

results_k0.25 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_0;25.csv")
results_k0.50 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_0;50.csv")
results_k0.75 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_0;75.csv")
results_k1.00 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_1;00.csv")
results_k1.50 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_1;50.csv")
results_k2.00 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_2;00.csv")


results_h1_k0.25 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_0;25.csv")
results_h1_k0.50 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_0;50.csv")
results_h1_k0.75 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_0;75.csv")
results_h1_k1.00 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_1;00.csv")
results_h1_k1.50 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_1;50.csv")
results_h1_k2.00 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_2;00.csv")










h0_allocation = data.frame(rbind((results_k0.25[1,]),
                                 (results_k0.50[1,]),
                                 (results_k0.75[1,]),
                                 (results_k1.00[1,]),
                                 (results_k1.50[1,]),
                                 (results_k2.00[1,])
                                 ))

h0_allocation = cbind("k" = c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00), h0_allocation)
write.csv(h0_allocation, "Results/allocation_ratio_h0_summary.csv", row.names = FALSE)






h1_allocation = data.frame(rbind((results_h1_k0.25[1,]),
                                 (results_h1_k0.50[1,]),
                                 (results_h1_k0.75[1,]),
                                 (results_h1_k1.00[1,]),
                                 (results_h1_k1.50[1,]),
                                 (results_h1_k2.00[1,])
                                  ))

h1_allocation = cbind("k" = c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00), h1_allocation)
write.csv(h1_allocation, "Results/allocation_ratio_h1_summary.csv", row.names = FALSE)







# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)
footnote = "Allocation Ratio: Blinded | H0 Holds | Alpha Estimate | sigmaTRUE = 6.5 | sigmaAssumed = 6.0"
plot_conf_intrv(h0_allocation$Alpha.Re.estimate, 
                h0_allocation$Confidence.Interval.High, 
                h0_allocation$Confidence.Interval.Low,
                footnote,
                0)




# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)
footnote = "Allocation Ratio: Blinded | H1 Holds | Beta Estimate | sigmaTRUE = 6.5 | sigmaAssumed = 6.0"
plot_conf_intrv(h1_allocation$Beta.Re.estimate, 
                h1_allocation$Confidence.Interval.Beta.High, 
                h1_allocation$Confidence.Interval.Beta.Low,
                footnote,
                1)



# latex table
cols = c(1, 5, 6, 7, 8)
print(xtable(h0_allocation[,cols], align = c("|c","|c","|c","|c","|c","|c|"), digits=c(0,2,3,0,0,5)), include.rownames = FALSE)



# latex table
cols = c(1, 5, 6, 7, 8)
print(xtable(h1_allocation[,cols], align = c("|c","|c","|c","|c","|c","|c|"), digits=c(0,2,3,0,0,5)), include.rownames = FALSE)





## APPENDIX ----


### Method 1: Base Method ----


#### H0 holds | blinded ----

met1_blinded_h0_v2 = read.csv("Results/met1_blinded_h0.csv")
met1_blinded_h0_distribution_v2 = read.csv("Results/met1_blinded_h0_distribution.csv")

# Histogram
ggplot(data=met1_blinded_h0_distribution_v2, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met1_blinded_h0_distribution_v2$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 1: H0 holds", subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")



#### H1 holds | blinded ----


met1_blinded_h1_v2 = read.csv("Results/met1_blinded_h1.csv")
met1_blinded_h1_distribution_v2 = read.csv("Results/met1_blinded_h1_distribution.csv")



# Histogram
ggplot(data=met1_blinded_h1_distribution_v2, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met1_blinded_h1_distribution_v2$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 1: H1 holds", subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")


#### H1 holds | partially unblinded (control) ----


met1_part_unblinded_C_h1_v2 = read.csv("Results/met1_part_unblinded_C_h1.csv")
met1_part_unblinded_C_h1_distribution_v2 = read.csv("Results/met1_part_unblinded_C_h1_distribution.csv")



# Histogram
ggplot(data=met1_part_unblinded_C_h1_distribution_v2, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met1_part_unblinded_C_h1_distribution_v2$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 1: H1 holds (Partially Unblinded - Control)", 
       subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")




#### H1 holds | partially unblinded (treatment) ----

met1_part_unblinded_T_h1_v2 = read.csv("Results/met1_part_unblinded_T_h1.csv")
met1_part_unblinded_T_h1_distribution_v2 = read.csv("Results/met1_part_unblinded_T_h1_distribution.csv")


# Histogram
ggplot(data=met1_part_unblinded_T_h1_distribution_v2, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met1_part_unblinded_T_h1_distribution_v2$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 1: H1 holds (Partially Unblinded - Treatment)", 
       subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")






### Method 2: Betensky Method ----



#### H0 holds | blinded ----

met2_blinded_h0_v3 = read.csv("Results/met2/met2_blinded_h0_v3.csv")
met2_blinded_h0_distribution_v3 = read.csv("Results/met2/met2_blinded_h0_distribution_v3.csv")


# Histogram
ggplot(data=met2_blinded_h0_distribution_v3, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met2_blinded_h0_distribution_v3$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 2: H0 holds", subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")


#### H1 holds | blinded ----


met2_blinded_h1_v3 = read.csv("Results/met2/met2_blinded_h1_v3.csv")
met2_blinded_h1_distribution_v3 = read.csv("Results/met2/met2_blinded_h1_distribution_v3.csv")


# Histogram
ggplot(data=met2_blinded_h1_distribution_v3, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met2_blinded_h1_distribution_v3$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 2: H1 holds (Blinded)", subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")





### Method 3: Bristol Method ----



#### H0 holds | blinded ----

met3_blinded_h0_v2 = read.csv("Results/met3_blinded_h0.csv")
met3_blinded_h0_distribution_v2 = read.csv("Results/met3_blinded_h0_distribution.csv")



# Histogram
ggplot(data=met3_blinded_h0_distribution_v2, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met3_blinded_h0_distribution_v2$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 3: H0 holds", subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")




#### H1 holds | blinded ----


met3_blinded_h1_v2 = read.csv("Results/met3_blinded_h1.csv")
met3_blinded_h1_distribution_v2 = read.csv("Results/met3_blinded_h1_distribution.csv")


# Histogram
ggplot(data=met3_blinded_h1_distribution_v2, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met3_blinded_h1_distribution_v2$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 3: H1 holds (Blinded)", subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")








### Method 4: Kieser Method ----



#### H0 holds | blinded ----


met4_blinded_h0_v2 = read.csv("Results/met4_blinded_h0.csv")
met4_blinded_h0_distribution_v2 = read.csv("Results/met4_blinded_h0_distribution.csv")


# Histogram
ggplot(data=met4_blinded_h0_distribution_v2, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met4_blinded_h0_distribution_v2$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 4: H0 holds", subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")



#### H1 holds | blinded ----


met4_blinded_h1_v2 = read.csv("Results/met4_blinded_h1.csv")
met4_blinded_h1_distribution_v2 = read.csv("Results/met4_blinded_h1_distribution.csv")


# Histogram
ggplot(data=met4_blinded_h1_distribution_v2, aes(Combination_6)) + 
  geom_histogram(aes(y =..density..)) +
  geom_density(col="blue", lwd = 2) +
  geom_vline(xintercept =  mean(met4_blinded_h1_distribution_v2$Combination_6), linetype="dashed", 
             color = "red", size=1) +
  labs(title = "Re-estimated Sample Size | Method 4: H1 holds (Blinded)", subtitle = "sigmaTRUE = 5.0 | sigma = 4.5")+
  xlab("Re-estimated Sample Size")




### ALLOCATION RATIO ----

results_k0.25 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_0;25.csv")
results_k0.50 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_0;50.csv")
results_k0.75 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_0;75.csv")
results_k1.00 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_1;00.csv")
results_k1.50 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_1;50.csv")
results_k2.00 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h0_k_2;00.csv")


results_h1_k0.25 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_0;25.csv")
results_h1_k0.50 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_0;50.csv")
results_h1_k0.75 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_0;75.csv")
results_h1_k1.00 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_1;00.csv")
results_h1_k1.50 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_1;50.csv")
results_h1_k2.00 = read.csv("Results/v2_allocation_ratio/allocation_study_blinded_h1_k_2;00.csv")










h0_allocation = data.frame(rbind((results_k0.25[5,]),
                                 (results_k0.50[5,]),
                                 (results_k0.75[5,]),
                                 (results_k1.00[5,]),
                                 (results_k1.50[5,]),
                                 (results_k2.00[5,])
))

h0_allocation = cbind("k" = c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00), h0_allocation)






h1_allocation = data.frame(rbind((results_h1_k0.25[5,]),
                                 (results_h1_k0.50[5,]),
                                 (results_h1_k0.75[5,]),
                                 (results_h1_k1.00[5,]),
                                 (results_h1_k1.50[5,]),
                                 (results_h1_k2.00[5,])
))

h1_allocation = cbind("k" = c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00), h1_allocation)







# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)
footnote = "Allocation Ratio: Blinded | H0 Holds | Alpha Estimate | sigmaTRUE = 5.0 | sigmaAssumed = 4.5"
plot_conf_intrv(h0_allocation$Alpha.Re.estimate, 
                h0_allocation$Confidence.Interval.High, 
                h0_allocation$Confidence.Interval.Low,
                footnote,
                0)




# Confidence Interval plot. plot_conf_intrv(estimate, CI_high, CI_low)
footnote = "Allocation Ratio: Blinded | H1 Holds | Beta Estimate | sigmaTRUE = 5.0 | sigmaAssumed = 4.5"
plot_conf_intrv(h1_allocation$Beta.Re.estimate, 
                h1_allocation$Confidence.Interval.Beta.High, 
                h1_allocation$Confidence.Interval.Beta.Low,
                footnote,
                1)



# latex table
cols = c(1, 5, 6, 7, 8)
print(xtable(h0_allocation[,cols], align = c("|c","|c","|c","|c","|c","|c|"), digits=c(0,2,3,0,0,5)), include.rownames = FALSE)



# latex table
cols = c(1, 5, 6, 7, 8)
print(xtable(h1_allocation[,cols], align = c("|c","|c","|c","|c","|c","|c|"), digits=c(0,2,3,0,0,5)), include.rownames = FALSE)













