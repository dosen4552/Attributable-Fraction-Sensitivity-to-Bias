a1 = 802
b1 = 166
c1 = 699
d1 = 269
a2 = 319
b2 = 180
c2 = 303
d2 = 196

attributable_effect = function(a, b, c, d , Gamma = 1, Theta = 1){
  A = a
  r = a + b
  m = a + c
  n = a + b + c + d
  f = function(x) x^2 * (Gamma * Theta - 1) - 
    x * ( (Gamma * Theta - 1 ) * (m + r) + n ) + Gamma * Theta * r * m
  E = uniroot(f, lower = max(0,r+m-n), upper = min(r,m))$root
  V = 1 / (1/E + 1/(r-E) + 1/(m-E) + 1/(n-r-m+E))
  return( 1 - pnorm((A - E)/sqrt(V)) )
}


a.vec=seq(0,a1+a2,1)
pval.right.vec = rep(0,length(a.vec))
pval.left.vec = rep(0,length(a.vec))
pval.vec=rep(0,length(a.vec))
gamma = 1.3
theta = 1.15

for(i in 1:length(a.vec)){
  a=a.vec[i]
  pval.vec[i]=attributable_effect(a1+a2-a,b1+b2,c1+c2+a,d1+d2, 
                                  Gamma = gamma, Theta = theta)
}
lci.total=min(a.vec[which(pval.vec>=.05)])


# Attributable effect using two tables
library(sensitivitymv)
a.clear.vec=seq(0,a1,1)
pval.clear.vec=rep(0,length(a.clear.vec))
for(i in 1:length(a.clear.vec)){
  a=a.clear.vec[i]
  pval.clear.vec[i]=attributable_effect(a1-a,b1,c1+a,d1, 
                                        Gamma = gamma, Theta = theta)
}

a.unclear.vec=seq(0,a2,1)
pval.unclear.vec=rep(0,length(a.unclear.vec))
for(i in 1:length(a.unclear.vec)){
  a=a.unclear.vec[i]
  pval.unclear.vec[i]=attributable_effect(a2-a,b2,c2+a,d2, 
                                          Gamma = gamma, Theta = theta)
}

acomb.vec=rep(0,length(a.clear.vec)*length(a.unclear.vec))
pval.comb.vec=rep(0,length(a.clear.vec)*length(a.unclear.vec))
pval.truncated=rep(0,length(a.clear.vec)*length(a.unclear.vec))
pval.bonferroni=rep(0,length(a.clear.vec)*length(a.unclear.vec))
count=1
for(i in 1:length(a.clear.vec)){
  for(j in 1:length(a.unclear.vec)){
    acomb.vec[count]=a.clear.vec[i]+a.unclear.vec[j]
    pval.comb.vec[count]=1-pchisq(-2*(log(pval.clear.vec[i]) + log(pval.unclear.vec[j]) )
                                  , 4, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    pval.truncated[count] = truncatedP(c(pval.clear.vec[i],pval.unclear.vec[j]),trunc=0.20)
    pval.bonferroni[count] = 2*min(c(pval.clear.vec[i],pval.unclear.vec[j]))
    count=count+1
  }
}

lci.comb=min(acomb.vec[which(pval.comb.vec>=.05)])
lci.truncated=min(acomb.vec[which(pval.truncated>=.05)])
lci.bonferroni=min(acomb.vec[which(pval.bonferroni>=.05)])


print(lci.total/1121)
print(lci.comb/1121)
print(lci.truncated/1121)
print(lci.bonferroni/1121)



print(lci.total)
print(lci.comb)
print(lci.truncated)
print(lci.bonferroni)


