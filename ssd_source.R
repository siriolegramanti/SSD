# extractorRData extracts an object from a .RData file created with save()
# Input: .RData filename and object name, both in " "

extractorRData <- function(file, object) {
  E <- new.env()
  load(file=file, envir=E)
  return(get(object, envir=E, inherits=F))
}


# F2 is needed to compute the test statistics in Barrett and Donald (2003) 
# for  2nd order stochastic dominance;
# Input:  a grid x of length N_grid and a sample a of length n,
# Output: a vector of length N_grid

F2 = function(x,a){
  order=2
  k = 1/gamma(order)
  M = matrix(x,length(x),1)%*%matrix(1,1,length(a))-matrix(1,length(x),1)%*%matrix(a,1,length(a))
  return(k*rowMeans(pmax(M,0))^(order-1))
}


# ssd_exp implements 3 tests (sup-based, integral-based, Barrett and Donald's)
# and saves both settings and results in the specified path
# Inputs:
# res_path = path where results will be saved
# a_path = path where data from the 1st distribution are stored 
# b_path = path where data from the 2nd distribution are stored
# n = number of data points to be taken from a_path and b_path
# (both data are stored as a matrix of size N_exp x max_n, 
# with max_n being the largest possible sample size)
# alpha = level of the test
# dependence = T if the two data distributions are dependent
# test_sup = T if the sup-based test is to be performed
# test_int = T if the integral-based test is to be performed
# test_Bar = T if Barrett and Donald's test is to be performed
# N_exp = number of experiments
# N_boot = number of bootstrap replicates
# N_grid = number of grid points for Barrett's test
# c = data shift

ssd_exp = function(res_path,a_path,b_path,n,alpha,
                   dependence=F,test_sup=T,test_int=T,test_Bar=T,
                   N_exp=500,N_boot=500,N_grid=100,c=10^(-4)){

  #load data
  A = extractorRData(a_path,"my_data") + c
  B = extractorRData(b_path,"my_data") + c
  if (n>ncol(A)|n>ncol(B))
    stop("n must be smaller than both ncol(A) and ncol(B)")

  # vector/matrices to store test statistics
  T_sup = rep(NA,N_exp)
  T_int = rep(NA,N_exp)
  T_Barrett = rep(NA,N_exp)
  Ts_sup = matrix(NA,N_exp,N_boot)
  Ts_int = matrix(NA,N_exp,N_boot)
  Ts_Barrett = matrix(NA,N_exp,N_boot)
  
  # useful grids
  xa = c(1:n)/n
  xm = xa-1/(2*n)
  
  # EXPERIMENTS
  st = Sys.time()
  for (i in 1:N_exp){
    print(i)
    a = A[i,1:n]
    b = B[i,1:n]
    # get the Lorenz curve of a
    ya = c(0,cumsum(sort(a))/n)
    Lor_a = stepfun(xa,ya)
    # get the INVERSE Lorenz curve of b
    xb = cumsum(sort(b))/n
    yb = c(0:n)/n
    if (mean(a)>mean(b)){
      xb = c(xb,mean(a))
      yb = c(yb,1)
    }
    Inv_Lor_b = stepfun(xb,yb)
    Phi = function(t){Inv_Lor_b(Lor_a(t))}
    # sup test statistics
    if (test_sup)
      T_sup[i] = sqrt(n)*max(0-Phi(0),xa-Phi(xa))
    # integral test statistics
    if (test_int){
      h = function(t){pmax(0,t-Phi(t))}
      T_int[i]=sqrt(n)*mean(h(xm))
    }
    # Barrett and Donald's test statistics
    if (test_Bar){
      grid = seq(from=min(c(a,b)),to=max(c(a,b)),length.out=N_grid)
      T_Barrett[i] = sqrt(n)*max(F2(grid,a)-F2(grid,b))
    }
    # BOOTSTRAP REPLICATES
    for (j in 1:N_boot){
      if (!dependence){
        # sample (with replacement) as from a and bs from b INDEPENDENTLY
        as = sample(a,size=n,replace=TRUE) 
        bs = sample(b,size=n,replace=TRUE) 
      }else{
        # sample (with replacement) n couples (a[i],b[i]) from (a[1:n],b[1:n])
        is = sample(c(1:n),size=n,replace=TRUE) #indexes
        as = a[is]
        bs = b[is]
      }
      # get the Lorenz curve of as
      yas = c(0,cumsum(sort(as))/n)
      Lor_as = stepfun(xa,yas)
      # get the INVERSE Lorenz curve of bs
      xbs = cumsum(sort(bs))/n
      ybs = c(0:n)/n
      if (mean(as)>mean(bs)){
        xbs = c(xbs,mean(as))
        ybs = c(ybs,1)
      }
      Inv_Lor_bs = stepfun(xbs,ybs)
      Phis = function(t){Inv_Lor_bs(Lor_as(t))}
      # sup test statistics
      if (test_sup)
        Ts_sup[i,j] = sqrt(n)*max(Phi(0)-Phis(0),Phi(xa)-Phis(xa))
      # integral test statistics
      if (test_int){
        hs = function(t){pmax(0,Phi(t)-Phis(t))}
        Ts_int[i,j] = sqrt(n)*mean(hs(xm))
      }
      # Barrett and Donald's test statistics
      if (test_Bar){
        grid = seq(from=min(c(as,bs)),to=max(c(as,bs)),length.out=N_grid)
        Ts_Barrett[i,j] = sqrt(n)*max(F2(grid,a)-F2(grid,b)-F2(grid,as)+F2(grid,bs))
      }
    }
  }
  # running time
  rt = Sys.time()-st
  
  # p-values
  p_int = rep(NA,N_exp)
  p_sup = rep(NA,N_exp)
  p_Barrett = rep(NA,N_exp)
  for (i in 1:N_exp){
    p_int[i]=mean(Ts_int[i,]>T_int[i])
    p_sup[i]=mean(Ts_sup[i,]>T_sup[i])
    p_Barrett[i]=mean(Ts_Barrett[i,]>T_Barrett[i])
  }

  # power of the test
  pow_int=mean(p_int<alpha)
  pow_sup=mean(p_sup<alpha)
  pow_Barrett=mean(p_Barrett<alpha)
  
  save(a_path,b_path,n,alpha,dependence,N_exp,N_boot,N_grid,c,
       rt,T_sup,Ts_sup,T_int,Ts_int,T_Barrett,Ts_Barrett,
       p_int,p_sup,p_Barrett,pow_int,pow_sup,pow_Barrett,
       file=res_path)
}