# Project with Kan
# Attributable effect for total table (cohort studies)
a1 = 802
b1 = 166
c1 = 699
d1 = 269
a2 = 319
b2 = 180
c2 = 303
d2 = 196

library('rbounds')
library(metap)
a.vec=seq(0,a1+a2,1)
pval.right.vec = rep(0,length(a.vec))
pval.left.vec = rep(0,length(a.vec))
pval.vec=rep(0,length(a.vec))
gamma = 1


for(i in 1:length(a.vec)){
  temp.table=matrix(c( a1+a2,b1+b2,c1+c2,d1+d2),ncol=2)
  a=a.vec[i]
  adjusted.table=temp.table+matrix(c(-a,0,a,0),ncol=2)
  pval.vec[i]=fisher.test(adjusted.table)$p.val
  
  #pval.right.vec[i]=FisherSens(totalN = sum(adjusted.table), 
                         #treatedN = sum(adjusted.table[1,]), 
                         #totalSuccesses = sum(adjusted.table[,1]), 
                         #treatedSuccesses = sum(adjusted.table[1,1]), Gammas = gamma)[1,3]
  #pval.left.vec[i] = 1 - FisherSens(totalN = sum(adjusted.table), 
                                #treatedN = sum(adjusted.table[1,]), 
                                #totalSuccesses = sum(adjusted.table[,1]), 
                                #treatedSuccesses = sum(adjusted.table[1,1]), Gammas = gamma)[1,2]
  #pval.vec[i] = 2 * min(c(pval.right.vec[i], pval.left.vec[i]))
}
lci.total=min(a.vec[which(pval.vec>=.05)])
uci.total=max(a.vec[which(pval.vec>=.05)])


# Attributable effect using two tables
library(sensitivitymv)
a.clear.vec=seq(0,a1,1)
pval.clear.vec=rep(0,length(a.clear.vec))
#pval.clear.right.vec = rep(0,length(a.clear.vec))
#pval.clear.left.vec = rep(0,length(a.clear.vec))
for(i in 1:length(a.clear.vec)){
  temp.table=matrix(c(a1,b1,c1,d1),ncol=2)
  a=a.clear.vec[i]
  adjusted.table=temp.table+matrix(c(-a,0,a,0),ncol=2)
  pval.clear.vec[i]=fisher.test(adjusted.table)$p.val
  #pval.clear.right.vec[i]=FisherSens(totalN = sum(adjusted.table), 
                               #treatedN = sum(adjusted.table[1,]), 
                               #totalSuccesses = sum(adjusted.table[,1]), 
                               #treatedSuccesses = sum(adjusted.table[1,1]), Gammas = gamma)[1,3]
  #pval.clear.left.vec[i] = 1 - FisherSens(totalN = sum(adjusted.table), 
                                    #treatedN = sum(adjusted.table[1,]), 
                                    #totalSuccesses = sum(adjusted.table[,1]), 
                                    #treatedSuccesses = sum(adjusted.table[1,1]), Gammas = gamma)[1,2]
  #pval.clear.vec[i] = 2 * min(c(pval.right.vec[i], pval.left.vec[i]))
}

a.unclear.vec=seq(0,a2,1)
pval.unclear.vec=rep(0,length(a.unclear.vec))
#pval.unclear.right.vec=rep(0,length(a.unclear.vec))
#pval.unclear.left.vec=rep(0,length(a.unclear.vec))
for(i in 1:length(a.unclear.vec)){
  temp.table=matrix(c(a2,b2,c2,d2),ncol=2)
  a=a.unclear.vec[i]
  adjusted.table=temp.table+matrix(c(-a,0,a,0),ncol=2)
  pval.unclear.vec[i]=fisher.test(adjusted.table)$p.val
  #pval.unclear.right.vec[i]=FisherSens(totalN = sum(adjusted.table), 
                                     #treatedN = sum(adjusted.table[1,]), 
                                     #totalSuccesses = sum(adjusted.table[,1]), 
                                     #treatedSuccesses = sum(adjusted.table[1,1]), Gammas = gamma)[1,3]
  #pval.unclear.left.vec[i] = 1 - FisherSens(totalN = sum(adjusted.table), 
                                          #treatedN = sum(adjusted.table[1,]), 
                                          #totalSuccesses = sum(adjusted.table[,1]), 
                                          #treatedSuccesses = sum(adjusted.table[1,1]), Gammas = gamma)[1,2]
  #pval.unclear.vec[i] = 2 * min(c(pval.right.vec[i], pval.left.vec[i]))
}

acomb.vec=rep(0,length(a.clear.vec)*length(a.unclear.vec))
pval.comb.vec=rep(0,length(a.clear.vec)*length(a.unclear.vec))
count=1
for(i in 1:length(a.clear.vec)){
  for(j in 1:length(a.unclear.vec)){
    acomb.vec[count]=a.clear.vec[i]+a.unclear.vec[j]
    #pval.comb.vec[count]=truncatedP(c(pval.clear.vec[i],pval.unclear.vec[j]),trunc=0.05)
    pval.comb.vec[count]=1-pchisq(-2*(log(pval.clear.vec[i]) + log(pval.unclear.vec[j]) )
                                  , 4, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    #alpha = 0.0255
    #w = (pval.clear.vec[i]^(pval.clear.vec[i] < alpha)) * (pval.unclear.vec[j]^(pval.unclear.vec[j] < alpha))
    #pval.comb.vec[count]=sum(dbinom(1:2,2,alpha)*(1-pgamma(-log(w/(alpha^(1:2))),1:2)))
    #pval.comb.vec[count]=2*min(c(pval.clear.vec[i],pval.unclear.vec[j]))
    #pval.comb.vec[count]=sumz(c(pval.clear.vec[i],pval.unclear.vec[j]),c(1163/2353,1190/2353))$p
    count=count+1
  }
}

lci.two=min(acomb.vec[which(pval.comb.vec>=.05)])
uci.two=max(acomb.vec[which(pval.comb.vec>=.05)])


print(lci.total)
print(uci.total)
print(lci.two)
print(uci.two)





