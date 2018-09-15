appGPEV <-
  function(train_data,test_data,delta,J,L,theta) { 
    
    print("method = appGPEV")
    ## start the time 
    ptm <- proc.time()   
    
    
    sigma=.2  ## sd
    
    loglik.fun=function(w,x,u,a) sum(dnorm(z,sqrt(2)*colSums(a*cos(outer(w,x)+u))/sqrt(J),sd=sigma,log=TRUE))
    
    
    # mu0=theta[1]
    # alpha0=theta[2]
    # ka0=theta[3] 
    # atau=theta[4]  # hyperparameter of \tau 
    # btau=theta[5]  # hyperparameter of \tau
    
    mu0=0
    alpha0=1
    ka0=1 
    atau=1  # hyperparameter of \tau 
    btau=1  # hyperparameter of \tau
    
    a0=theta[1]
    b0=theta[2] # lambda \sim gamma(a0,b0)
    lambda0=theta[3]
    K=25
    s=theta[4]
    
    n=ncol(train_data)
    locs=as.numeric(train_data[1,])
    locs.err=as.numeric(train_data[2,])
    z=as.numeric(train_data[3,])
    zt=as.numeric(train_data[4,])
    n.p=ncol(test_data)
    predlocs=as.numeric(test_data[1,])
    z1=as.numeric(test_data[2,])
    
    w=u=a=matrix(nrow=J,ncol=L)
    x=matrix(nrow=n,ncol=L)
    BPI=matrix(nrow=L,ncol=K) # pi--stick-breaking 
    BZ=matrix(nrow=L,ncol=n) # cluster labels
    mu=matrix(nrow=L,ncol=K) 
    tau=matrix(nrow=L,ncol=K)
    t=matrix(nrow=n,ncol=K)
    lambda=c();lambda[1]=lambda0
    b=c(); b[1]=b0
    V0=c(rbeta(K-1,1,alpha0), 1)
    tau[1,]=rgamma(K,atau,btau)
    mu[1,]=rnorm(K,mu0,sd=sqrt(ka0/tau[1,]))
    
    for(h in 1:K-1){
      BPI[1,h]=V0[h]*prod(1-V0[1:h-1])
    }
    BPI[1,K]=prod(1-V0[1:K-1])
    
    
    
    w[,1]=w.cur=rep((1:(J/2))*2*pi/n,2)
    u[,1]=u.cur=c(rep(pi/2,J/2),rep(3*pi/2,J/2))
    a[,1]=a.cur=rnorm(J,0,1)
    x[,1]=x.cur=locs.err
    loglik.cur=loglik.fun(w.cur,x.cur,u.cur,a.cur)
    
    
    ## mcmc 
    for(l in 2:L) { # l=2
     #cat("iter=", l, "\n")
      for(j in 1:J) { # j=1
        w.prop=w.cur; w.prop[j]=1.5*abs(rnorm(1,0,sd=8*sqrt(1/lambda[l-1]))/s)
        loglik.prop=loglik.fun(w.prop,x.cur,u.cur,a.cur)
        if(loglik.prop-loglik.cur>log(runif(1))){
          w.cur=w.prop; loglik.cur=loglik.prop 
        }
        u.prop=u.cur; u.prop[j]=(rnorm(1,u.cur[j],sd=.5)) %% (2*pi)
        loglik.prop=loglik.fun(w.cur,x.cur,u.prop,a.cur)
        if(loglik.prop-loglik.cur>log(runif(1))){
          u.cur=u.prop; loglik.cur=loglik.prop
        } 
        
        
        w[,l]=w.cur; u[,l]=u.cur; 
      }
      
      # update a_j
      phi=sqrt(2)*cos(outer(w[,l],x[,l-1])+u[,l])/sqrt(J)/sigma
      r=rnorm(J,0,1); d=rnorm(n,0,1)
      v=t(phi)%*%r+d; m=solve(t(phi)%*%phi+diag(n),z/sigma-v)
      a[,l]=r+phi%*%m; a.cur=a[,l]
      
      ## DPMM updating f_X
      # update cluster labels Z_i
      for(j in 1:n){
        t[j,]=dnorm(x.cur[j],mu[l-1,],sd=1/sqrt(tau[l-1,]))
        pz=BPI[l-1,]*t[j,]/sum(BPI[l-1,]*t[j,])
        BZ[l-1,j]=sample(c(1:K),1,prob=pz)
      }
      
      # update \mu and \tau 
      yhat=c();sh=c();nh=c()
      for(h in 1:K){
        nh[h]=sum((BZ[l-1,]==h))
        if(nh[h]>0){
          yhat[h]=sum(x.cur[which(BZ[l-1,]==h)])/nh[h]
          sh[h]=sum((x.cur[which(BZ[l-1,]==h)]-yhat[h])^2)
        }
        else{
          yhat[h]=0; sh[h]=0 
        }
      }
      ka=1/(1/ka0+nh) # update kappa
      muh=ka*(mu0/ka0 + nh*yhat) 
      a_tau=atau+nh/2
      b_tau=btau + 0.5*(sh+nh*(yhat-mu0)^2/(1+ka0*nh))
      tau[l,]=rgamma(K,a_tau,b_tau)
      mu[l,]=rnorm(K,muh,sd=sqrt(ka/(tau[l,])))
      
      # update \pi_h, h=1,...,K
      N=1+nh;Nb=rcumsum(nh)+alpha0
      V=c(rbeta(K-1,N[1:K-1],Nb[-1]),1)
      for(h in 1:K-1){
        BPI[l,h]=V[h]*prod(1-V[1:h-1])
      }
      BPI[l,K]=prod(1-V[1:K-1])
      
      ## updating x
      tau.cur=tau[l,];mu.cur=mu[l,]
      for (j in 1:n){
        id=sample(c(1:K),1,prob=BPI[l,])
        x.prop=x.cur; x.var=1/delta+tau.cur[id]
        x.prop[j]=rnorm(1,(locs.err[j]/delta+mu.cur[id]*tau.cur[id])/x.var,sd=sqrt(1/x.var))
        loglik.prop=loglik.fun(w.cur,x.prop,u.cur,a.cur)
        if(loglik.prop-loglik.cur >log(runif(1))){
          x.cur=x.prop; loglik.cur=loglik.prop
        }
      }
      x[,l]=x.cur
      
      
      # update lambda
      b[l]=b0/(1+b0*sum(w.cur^2)/4)
      lambda[l]=rgamma(1,a0,scale = 1/b[l])
      
      
    }
    
    
    w.post=t(w[,seq(200,L,length=500)])
    u.post=t(u[,seq(200,L,length=500)])
    a.post=t(a[,seq(200,L,length=500)])
    x.post=rowMeans(x[,seq(200,L,length=500)])
    
    mcmc_w = w[,(L-200+1):L]
    mcmc_u = u[,(L-200+1):L]
    mcmc_a = a[,(L-200+1):L]
    mcmc_x = x[,(L-200+1):L]
    
   # mcmc_sample = list(mcmc_w, mcmc_u, mcmc_a, mcmc_x)
    
    # \hat{f}(\hat{x})
    bfs=array(dim=c(n,J,500)); for(j in 1:J){ for(l in 1:500) bfs[,j,l]=a.post[l,j]*cos(outer(w.post[l,j],x.post)+u.post[l,j])}
    y.post=sqrt(2/J)*apply(bfs,c(1,3),sum)
    ypost=apply(y.post,1,mean)
    mse.train=sum((ypost-zt)^2)/n
    
    bfs.mc=array(dim=c(n,J,L))
    for(j in 1:J){ for(l in 1:L) bfs.mc[,j,l]=a[j,l]*cos(outer(w[j,l],x[,l])+u[j,l])}
    f=sqrt(2/J)*apply(bfs.mc,c(1,3),sum)
    r= abs(sweep(f,2,zt))
    rmax=apply(r,2,max);rq=quantile(rmax,probs=c(0.025,0.975))
    rqlb=rq[1];rqub=rq[2]
    
    
    # \hat{f}(x_0)
    bfs.locs=array(dim=c(n,J,500)); for(j in 1:J){ for(l in 1:500) bfs.locs[,j,l]=a.post[l,j]*cos(outer(w.post[l,j],locs)+u.post[l,j])}
    y.true=sqrt(2/J)*apply(bfs.locs,c(1,3),sum)
    ytrue=apply(y.true,1,mean)
    
    # out of sample preds
    bfs1=array(dim=c(n.p,J,500)); for(j in 1:J){ for(l in 1:500) bfs1[,j,l]=a.post[l,j]*cos(outer(w.post[l,j],predlocs)+u.post[l,j])}
    y1.pred=sqrt(2/J)*apply(bfs1,c(1,3),sum)
    ypred=apply(y1.pred,1,mean)
    MSE=sum((ypred-z1)^2)/n.p
    
    bfs1.mc = array(dim = c(n.p,J,1000))
    for(j in 1:J){ for(l in 1:1000) bfs1.mc[,j,l] = a[j,l+500]*cos(outer(w[j,l+500],predlocs)+u[j,l+500])}
    f1 = sqrt(2/J)*apply(bfs1.mc,c(1,3),sum)
    r1 = abs(sweep(f1,1,z1,'-'))
    rmax1 = apply(r1,2,max)
    r.pred = as.numeric(quantile(rmax1,probs = 0.95))
    lb = ypred - r.pred
    ub = ypred + r.pred
   
    sCI_test = rbind(lb, ub)
    pCI_test = apply(f1, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    
    ## Stop the time
    
    ptm <- proc.time() - ptm
    
    ## Processing time
    
    (minutes <- ptm[[3]]/60)
    
  
    result=list(list())
    
    result=list(ypost, ytrue, ypred, mse.train, MSE, sCI_test, pCI_test, minutes, mcmc_w,mcmc_u,mcmc_a,mcmc_x)
    
    
    return(result)
    
    
    
  }
