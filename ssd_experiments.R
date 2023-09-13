# The following code allows to reproduce the results in Table 6a.
# Table 6b is obtained by switching samples a and b.
# Tables 7 and 8 are obtained by changing the Singh-Maddala parameters.
# The other tables can be obtained by changing the data-generating process and,
# in case of dependent samples, by setting "dependence=T" in 'ssd_exp()'.

source("ssd_source.R")
library(VGAM) # to sample from the Singh-Maddala distribution

# sample sizes
ns = c(50,100,200,500,1000)
# number of datasets for each sample size
N_exp = 500

set.seed(123)

# first sample (a ~ F)
my_data=matrix(rsinmad(n=N_exp*max(ns),shape1.a=1.5,shape3.q=1.8),
               N_exp,max(ns),byrow=T)
save(my_data,file=paste0("a.RData"))

# second sample (b ~ G)
my_data=matrix(rsinmad(n=N_exp*max(ns),shape1.a=1,shape3.q=1.8),
               N_exp,max(ns),byrow=T)
save(my_data,file=paste0("b.RData"))

# experiments
for (n in ns){
  res_path=paste0("n_",n,".RData")
  print(res_path)
  ssd_exp(res_path,a_path="a.RData",b_path="b.RData",n,alpha=0.1,
    dependence=F,test_sup=T,test_int=T,test_Bar=T,
    N_exp,N_boot=500,N_grid=100,c=10^(-4))
}

#results table
res_table = matrix(NA,length(ns),4)
colnames(res_table) = c("n","sup","integral","Barrett")
for (i in 1:length(ns)){
  load(paste0("n_",ns[i],".RData"))
  res_table[i,]=c(n,pow_sup,pow_int,pow_Barrett)
}
res_table