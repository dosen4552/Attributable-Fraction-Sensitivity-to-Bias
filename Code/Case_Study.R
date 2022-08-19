library(sensitivitymv)

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

get_lower = function(vec){
  a = c(vec[1]) 
  for(i in 2: (length(vec)-1) ){
    if(vec[i] < vec[i-1] & vec[i] < vec[i+1]){
      a = c(a,vec[i])
    }
  }
  return(min(a))
}


get_upper = function(vec){
  a = c(vec[1]) 
  for(i in 2: (length(vec)-1) ){
    if(vec[i] < vec[i-1] & vec[i] < vec[i+1]){
      a = c(a,vec[i])
    }
  }
  return(max(a))
}

## Load Matched Data
path_in <- 'C:/Users/chenk/Dropbox/kan_dylan/attributable_effect_case_study/datasets/data_needed'

matched_data1 = read.csv(paste0(path_in, "/matched_data1.csv"))
matched_data2 = read.csv(paste0(path_in, "/matched_data2.csv"))




case1 = get_abcd(matched_data = matched_data1)
case2 = get_abcd(matched_data = matched_data2)



a1 = case1[1] # 1 
b1 = case1[2] # 86
c1 = case1[3] # 43
d1 = case1[4] # 3024
a2 = case2[1] # 1
b2 = case2[2] # 15
c2 = case2[3] # 21
d2 = case2[4] # 855



# Lower bound
# Merge result

t = 1.1
g_list = seq(1.0,1.4,0.02)
iter = 1.0
result.table = matrix(0, length(g_list), 5)

for (g in g_list) {
  a.vec=seq(0,b1 + b2 + c1 + c2,1)
  pval.vec=rep(0,length(a.vec))
  for(i in 1:length(a.vec)){
    a=a.vec[i]
    pval.vec[i]=attributable_effect_p_value_lower(a1 + a2,b1 + b2,c1 
                                                  + c2,d1 + d2,a,Gamma = t * g)
  }
  lci.total= min(a.vec[which(pval.vec>=.05)])
  
  
  # Two table result
  a.vec.1 = seq(0,b1 + c1,1)
  pval.vec.1 = rep(0,length(a.vec.1))
  
  for(i in 1:length(a.vec.1)){
    a=a.vec.1[i]
    pval.vec.1[i]=attributable_effect_p_value_lower(a1,b1,c1,d1,a,Gamma = t * g)
  }
  
  
  a.vec.2 = seq(0,b2 + c2,1)
  pval.vec.2 = rep(0,length(a.vec.2))
  
  for(i in 1:length(a.vec.2)){
    a=a.vec.2[i]
    pval.vec.2[i]=attributable_effect_p_value_lower(a2,b2,c2,d2,a,Gamma = t * g)
  }
  
  acomb.vec=rep(0,length(a.vec.1)*length(a.vec.2))
  pval.comb.vec.fisher=rep(0,length(a.vec.1)*length(a.vec.2))
  pval.comb.vec.zscore=rep(0,length(a.vec.1)*length(a.vec.2))
  pval.comb.vec.truncated005=rep(0,length(a.vec.1)*length(a.vec.2))
  pval.comb.vec.bonferroni=rep(0,length(a.vec.1)*length(a.vec.2))
  count=1
  
  w1 = 0.84 
  w2 = 0.16 
  
  
  for(i in 1:length(a.vec.1)){
    for(j in 1:length(a.vec.2)){
      acomb.vec[count]=a.vec.1[i]+a.vec.2[j]
      pval.comb.vec.fisher[count]=1-pchisq(-2*(log(pval.vec.1[i]) + 
                                                 log(pval.vec.2[j]) )
                              , 4, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      pval.comb.vec.zscore[count]=1 - pnorm( (( w1 *qnorm(1 - pval.vec.1[i]) + 
                                                  w2 * qnorm(1 - 
                                          pval.vec.2[i]))/sqrt(w1^2 + w2^2)) )
      pval.comb.vec.truncated005[count]=truncatedP(c(pval.vec.1[i],pval.vec.2[j])
                                                   , trunc = 0.10 )
      pval.comb.vec.bonferroni[count]=2*min(c(pval.vec.1[i],pval.vec.2[j]))
      count=count+1
    }
  }
  
  lci.fisher=min(acomb.vec[which(pval.comb.vec.fisher>=.05)])
  lci.zscore=min(acomb.vec[which(pval.comb.vec.zscore>=.05)])
  lci.truncated005=min(acomb.vec[which(pval.comb.vec.truncated005>=.05)])
  lci.bonferroni=min(acomb.vec[which(pval.comb.vec.bonferroni>=.05)])
  
  
  
  
  # Upper bound
  # Merge result
  a.vec=seq(0,b1 + b2 + c1 + c2,1)
  pval.vec=rep(0,length(a.vec))
  
  for(i in 1:length(a.vec)){
    a=a.vec[i]
    pval.vec[i]=attributable_effect_p_value_upper(a1 + a2,b1 + b2,c1 + c2,d1 
                                                  + d2,a,Gamma = t * g)
  }
  uci.total= max(a.vec[which(pval.vec>=.95)])
  
  
  # Two table result
  a.vec.1 = seq(0,b1 + c1,1)
  pval.vec.1 = rep(0,length(a.vec.1))
  
  for(i in 1:length(a.vec.1)){
    a=a.vec.1[i]
    pval.vec.1[i]=attributable_effect_p_value_upper(a1,b1,c1,d1,a,Gamma = t * g)
  }
  
  
  a.vec.2 = seq(0,b2 + c2,1)
  pval.vec.2 = rep(0,length(a.vec.2))
  
  for(i in 1:length(a.vec.2)){
    a=a.vec.2[i]
    pval.vec.2[i]=attributable_effect_p_value_upper(a2,b2,c2,d2,a,Gamma = t * g)
  }
  
  acomb.vec=rep(0,length(a.vec.1)*length(a.vec.2))
  pval.comb.vec.fisher=rep(0,length(a.vec.1)*length(a.vec.2))
  pval.comb.vec.zscore=rep(0,length(a.vec.1)*length(a.vec.2))
  pval.comb.vec.truncated005=rep(0,length(a.vec.1)*length(a.vec.2))
  pval.comb.vec.bonferroni=rep(0,length(a.vec.1)*length(a.vec.2))
  count=1
  
  for(i in 1:length(a.vec.1)){
    for(j in 1:length(a.vec.2)){
      acomb.vec[count]=a.vec.1[i]+a.vec.2[j]
      pval.comb.vec.fisher[count]=1-pchisq(-2*(log(pval.vec.1[i]) 
                                               + log(pval.vec.2[j]) )
                                , 4, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      pval.comb.vec.zscore[count]=1 - pnorm( (( w1*qnorm(1 - pval.vec.1[i]) + 
                                w2*qnorm(1 - pval.vec.2[i]))/sqrt(w1^2+w2^2)) )
      pval.comb.vec.truncated005[count]=truncatedP(c(pval.vec.1[i],pval.vec.2[j])
                                                   , trunc = 0.10 )
      pval.comb.vec.bonferroni[count]=2*min(c(pval.vec.1[i],pval.vec.2[j]))
      count=count+1
    }
  }
  
  uci.fisher=max(acomb.vec[which(pval.comb.vec.fisher>=.95)])
  uci.zscore=max(acomb.vec[which(pval.comb.vec.zscore>=.95)])
  uci.truncated005=max(acomb.vec[which(pval.comb.vec.truncated005>=.95)])
  uci.bonferroni=max(acomb.vec[which(pval.comb.vec.bonferroni>=.95)])
  
  
  # Only display lower bound
  
  
  result = c(lci.total/(a1+a2+b1+b2), lci.zscore/(a1+a2+b1+b2), 
             lci.fisher/(a1+a2+b1+b2),lci.truncated005/(a1+a2+b1+b2),
             lci.bonferroni/(a1+a2+b1+b2))
  
  
  result.table[iter,] = result
  iter = iter + 1
  print(iter)
}

result.table = cbind(g_list, rep(t,length(g_list)), result.table * 100)
colnames(result.table) = c('Gamma','Theta','Merged','Stouffer','Fisher'
                           ,'alpha = 0.10','Bonferroni')


print(result.table)



















