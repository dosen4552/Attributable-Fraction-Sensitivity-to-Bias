library(mvtnorm)
library(mvnfast)
library(bigmatch)
library(optmatch)
library(mvtnorm)
library(match2C)
library(sensitivitymv)
library(rbounds)
library(metap)


generate_data <- function(prob = c(4/5,1/5),d = 10, n_t = 500, ratio = 1.5){
  n_c = n_t * ratio
  Z = c(rep(1, n_t), rep(0, n_c))
  mean_treated = rep(0, d)
  mean_control = rep(0, d)
  sigma = diag(d)
  X_treated = rmvnorm(n_t, mean = mean_treated, sigma = sigma)
  X_control = rmvnorm(n_c, mean = mean_control, sigma = sigma)
  X = rbind(X_treated, X_control)
  # Generate two potential outcomes for each unit
  n = n_t + n_c
  R_c = sample(0:1,n,replace = TRUE,prob = prob)
  R_t = sample(0:1,n,replace = TRUE,prob = 1-prob)
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
  propens = glm(Z ~ X, family = 'binomial')$fitted.value
  dist_list_pscore_maha_soft = create_list_from_scratch(Z, X, exact = NULL, 
                                                        p = propens, 
                                                        caliper_low = 0.2*sd(propens), 
                                                        method = 'maha', 
                                                        penalty = 1000) 
  matching_output = match_2C_list(Z, dataset_in_dataframe, 
                                  dist_list_pscore_maha_soft, 
                                  dist_list_2 = NULL, 
                                  controls = 1)
}


main_simulation <- function(prob1 = c(0.6,0.4), prob2 = c(0.6,0.4)){
  n = 240
  d = 10
  g = seq(1,5,0.5)
  n_sim = 200
  freq.table = 0
  
  # Use merged table
  
  for (iter in 1:n_sim){
    
    data1 = generate_data(prob = prob1,d = d, n_t = n/2, ratio = 1.5)
    matched_data1 = opt_match(data1)$matched_data_in_order
    dataset1 = matched_data1[! is.na(matched_data1$matched_set),]
    
    a1 = dim(dataset1[dataset1$Z == 1 & dataset1$R == 1,])[1]
    b1 = dim(dataset1[dataset1$Z == 1 & dataset1$R == 0,])[1]
    c1 = dim(dataset1[dataset1$Z == 0 & dataset1$R == 1,])[1]
    d1 = dim(dataset1[dataset1$Z == 0 & dataset1$R == 0,])[1]
    
    
    data2 = generate_data(prob = prob2,d = d, n_t = n/2, ratio = 1.5)
    matched_data2 = opt_match(data2)$matched_data_in_order
    dataset2 = matched_data2[!is.na(matched_data2$matched_set),]
    
    a2 = dim(dataset2[dataset2$Z == 1 & dataset2$R == 1,])[1]
    b2 = dim(dataset2[dataset2$Z == 1 & dataset2$R == 0,])[1]
    c2 = dim(dataset2[dataset2$Z == 0 & dataset2$R == 1,])[1]
    d2 = dim(dataset2[dataset2$Z == 0 & dataset2$R == 0,])[1]
    
    
    
    temp.table=matrix(c( a1 + a2,b1 + b2,c1 + c2,d1 + d2),ncol=2)
    pval.total=FisherSens(sum(temp.table), colSums(temp.table)[1], rowSums(temp.table)[1], 
                          temp.table[1,1],Gamma = g)[,3]
    
    
    
    # Use two tables
    
    temp.table1=matrix(c(a1,b1,c1,d1),ncol=2)
    temp.table2=matrix(c(a2,b2,c2,d2),ncol=2)
    
    pval1=FisherSens(sum(temp.table1), colSums(temp.table1)[1], rowSums(temp.table1)[1], 
                     temp.table1[1,1],Gamma = g)[,3]
    
    pval2=FisherSens(sum(temp.table2), colSums(temp.table2)[1], rowSums(temp.table2)[1], 
                     temp.table2[1,1],Gamma = g)[,3]
    
    pval.fisher = rep(0, length(pval1))
    pval.truncated005 = rep(0, length(pval1))
    pval.truncated010 = rep(0, length(pval1))
    pval.truncated015 = rep(0, length(pval1))
    pval.truncated020 = rep(0, length(pval1))
    pval.bonferroni = rep(0, length(pval1))
    
    for (i in 1:(length(pval1))){
      pval.fisher[i] = 1-pchisq(-2*(log(pval1[i]) + log(pval2[i]) )
                                , 4, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      pval.truncated005[i] = truncated(c(pval1[i],pval2[i]), ptrunc = 0.05 )
      pval.truncated010[i] = truncated(c(pval1[i],pval2[i]), ptrunc = 0.10 )
      pval.truncated015[i] = truncated(c(pval1[i],pval2[i]), ptrunc = 0.15 )
      pval.truncated020[i] = truncated(c(pval1[i],pval2[i]), ptrunc = 0.20 )
      pval.bonferroni[i] = 2*min(c(pval1[i],pval2[i]))
    }
    
    
    
    
    pval.table = cbind(pval.total,pval.fisher,pval.truncated005, pval.truncated010, 
                       pval.truncated015, pval.truncated020, pval.bonferroni)
    freq.table = freq.table + (pval.table < .05)
    print(iter)
  }
  
  
  
  power.table = cbind(g,freq.table/n_sim) 
  colnames(power.table) = c('Gamma','Merged','Fisher','alpha = 0.05','alpha = 0.10',
                            'alpha = 0.15','alpha = 0.20','Bonferroni')
  
  
  print(t(power.table))
}




# Start simulation
main_simulation(c(0.7,0.3),c(0.6,0.4)) #Equal Effect

main_simulation(c(0.65,0.35),c(0.6,0.4)) #Slightly unequal Effect

main_simulation(c(0.8,0.2),c(0.55,0.45)) #Unequal Effect

main_simulation(c(0.8,0.2),c(0.5,0.5)) #Effect only in one table






