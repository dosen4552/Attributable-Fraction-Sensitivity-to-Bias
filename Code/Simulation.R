library(mvtnorm)
library(mvnfast)
library(bigmatch)
library(optmatch)
library(mvtnorm)
library(match2C)
library(sensitivitymv)
library(rbounds)
library(poolr)


# Attributable Effect P-values calculation for lower bound and upper bound
attributable_effect_p_value_lower = function(a,b,c,d,A,Gamma = 1){
  p = Gamma/(Gamma + 1)
  z = ( (a+b - A) - (a + (b+c-A)*p))/sqrt((b+c-A)*p * (1-p))
  p_value = 1 - pnorm(z)
  return(p_value)
}



attributable_effect_p_value_upper = function(a,b,c,d,A,Gamma = 1){
  p = 1/(Gamma + 1)
  z = ( (a+b - A) - (a + (b+c-A)*p))/sqrt((b+c-A)*p * (1-p))
  p_value = 1 - pnorm(z)
  return(p_value)
}


# Data generating process and matching
generate_data <- function(prob = c(4/5,1/5),d = 1, n_t = 500, ratio = 1.5){
  n_c = n_t * ratio
  Z = c(rep(1, n_t), rep(0, n_c))
  mean_treated = rep(0, d)
  mean_control = rep(0, d + 0.2)
  sigma = diag(d)
  X_treated = rmvnorm(n_t, mean = mean_treated, sigma = sigma)
  X_control = rmvnorm(n_c, mean = mean_control, sigma = sigma)
  X = rbind(X_treated, X_control)
  # Generate two potential outcomes for each unit
  n = n_t + n_c
  # R_t >= R_c
  R_t = sample(0:1,n,replace = TRUE,prob = 1-prob)
  R_c = rep(0,n)
  R_c[which(R_t == 1)] = sample(0:1,sum(R_t),replace = TRUE,prob = prob)
  R = (Z == 0)*R_c + (Z == 1)*R_t + 0
  # return a list of two objects: X and Z
  return(list(X = X, Z = Z, R = R))
}

opt_match <- function(dataset){
  Z = dataset$Z
  X = dataset$X
  R = dataset$R
  dataset_in_dataframe = data.frame(X, Z, R)
  # Estimate pscore
  propens = glm(R ~ X, family = 'binomial')$fitted.value
  dist_list_pscore_maha_soft = create_list_from_scratch(R, X, exact = NULL, 
                                                        p = propens, 
                                                        caliper_low = 0.2*sd(propens), 
                                                        method = 'maha', 
                                                        penalty = 1000) 
  matching_output = match_2C_list(R, dataset_in_dataframe, 
                                  dist_list_pscore_maha_soft, 
                                  dist_list_2 = NULL, 
                                  controls = 1)
}


get_abcd <- function(matched_data){
  a = 0
  b = 0
  c = 0
  d = 0
  for (i in unique(matched_data$matched_set)){
    temp = matched_data[matched_data$matched_set == i,]
    if(temp[temp$R == 1,]$Z == 1 & temp[temp$R == 0,]$Z == 1 ){
      a = a + 1
    }
    if(temp[temp$R == 1,]$Z == 1 & temp[temp$R == 0,]$Z == 0 ){
      b = b + 1
    }
    if(temp[temp$R == 1,]$Z == 0 & temp[temp$R == 0,]$Z == 1 ){
      c = c + 1
    }
    if(temp[temp$R == 1,]$Z == 0 & temp[temp$R == 0,]$Z == 0 ){
      d = d + 1
    }
  }
  return(c(a,b,c,d))
}


main_simulation <- function(prob1 = c(0.6,0.4), prob2 = c(0.6,0.4),theta = 1.0){
  n = 400
  d = 10
  g = seq(1,5,0.5)
  n_sim = 200
  freq.table = 0
  
  # Use merged table
  
  for (iter in 1:n_sim){
    
    data1 = generate_data(prob = prob1,d = d, n_t = n/2, ratio = 1.5)
    matched_data1 = opt_match(data1)$matched_data_in_order
    dataset1 = matched_data1[! is.na(matched_data1$matched_set),]
    

    a1 = get_abcd(matched_data = dataset1)[1]
    b1 = get_abcd(matched_data = dataset1)[2]
    c1 = get_abcd(matched_data = dataset1)[3]
    d1 = get_abcd(matched_data = dataset1)[4]
    
    data2 = generate_data(prob = prob2,d = d, n_t = n/2, ratio = 1.5)
    matched_data2 = opt_match(data2)$matched_data_in_order
    dataset2 = matched_data2[!is.na(matched_data2$matched_set),]
    
    a2 = get_abcd(matched_data = dataset2)[1]
    b2 = get_abcd(matched_data = dataset2)[2]
    c2 = get_abcd(matched_data = dataset2)[3]
    d2 = get_abcd(matched_data = dataset2)[4]
    
    
    
    pval.total=attributable_effect_p_value_lower(a1+a2, b1+b2, c1+c2,d1+d2, 0
                                                   , Gamma = g * theta)
    
    
    
    # Use two tables
    
    pval1=attributable_effect_p_value_lower(a1, b1, c1, d1, 0, Gamma = g * theta)
    
    pval2=attributable_effect_p_value_lower(a2, b2, c2, d2, 0, Gamma = g * theta)
    
    pval1[pval1 == 0] = 2.225074e-308
    pval2[pval2 == 0] = 2.225074e-308
    
    pval.fisher = rep(0, length(pval1))
    pval.zscore = rep(0, length(pval1))
    pval.truncated005 = rep(0, length(pval1))
    pval.truncated010 = rep(0, length(pval1))
    pval.truncated015 = rep(0, length(pval1))
    pval.truncated020 = rep(0, length(pval1))
    pval.bonferroni = rep(0, length(pval1))
    
    w1 = 0.6
    w2 = 0.4
    
    for (i in 1:(length(pval1))){
      pval.fisher[i] = 1-pchisq(-2*(log(pval1[i]) + log(pval2[i]) )
                                , 4, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      pval.zscore[i] = 1 - pnorm( ((w1 * qnorm(1 - pval1[i]) + 
                                  w2 * qnorm(1 - pval2[i]))/sqrt(w1^2 + w2^2)) )
      pval.truncated005[i] = truncatedP(c(pval1[i],pval2[i]), trunc = 0.05 )
      pval.truncated010[i] = truncatedP(c(pval1[i],pval2[i]), trunc = 0.10 )
      pval.truncated015[i] = truncatedP(c(pval1[i],pval2[i]), trunc = 0.15 )
      pval.truncated020[i] = truncatedP(c(pval1[i],pval2[i]), trunc = 0.20 )
      pval.bonferroni[i] = 2*min(c(pval1[i],pval2[i]))
    }
    
    
    
    
    pval.table = cbind(pval.total,pval.zscore,pval.fisher,pval.truncated005, 
                       pval.truncated010, 
                       pval.truncated015, pval.truncated020, pval.bonferroni)
    freq.table = freq.table + (pval.table < .05)
    print(iter)
  }
  
  
  
  power.table = cbind(g,freq.table/n_sim) 
  colnames(power.table) = c('Gamma','Merged','Z-score','Fisher','alpha = 0.05','alpha = 0.10',
                            'alpha = 0.15','alpha = 0.20','Bonferroni')
  
  print(t(power.table))
}



set.seed(1)
# Start simulation
# Theta = 1.0
main_simulation(c(0.7,0.5),c(0.6,0.4),theta = 1.0) #Equal Effect

main_simulation(c(0.65,0.35),c(0.6,0.4),theta = 1.0) #Slightly unequal Effect

main_simulation(c(0.8,0.2),c(0.55,0.45),theta = 1.0) #Unequal Effect

main_simulation(c(0.8,0.2),c(0.5,0.5),theta = 1.0) #Effect only in one table



# Theta = 1.1
main_simulation(c(0.7,0.5),c(0.6,0.4),theta = 1.1) #Equal Effect

main_simulation(c(0.65,0.35),c(0.6,0.4),theta = 1.1) #Slightly unequal Effect

main_simulation(c(0.8,0.2),c(0.55,0.45),theta = 1.1) #Unequal Effect

main_simulation(c(0.8,0.2),c(0.5,0.5),theta = 1.1) #Effect only in one table



# Theta = 1.2
main_simulation(c(0.7,0.5),c(0.6,0.4),theta = 1.2) #Equal Effect

main_simulation(c(0.65,0.35),c(0.6,0.4),theta = 1.2) #Slightly unequal Effect

main_simulation(c(0.8,0.2),c(0.55,0.45),theta = 1.2) #Unequal Effect

main_simulation(c(0.8,0.2),c(0.5,0.5),theta = 1.2) #Effect only in one table

