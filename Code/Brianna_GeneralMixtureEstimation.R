require(dplyr)
require(plyr)
require(nnet)

gen=function(Time,X,id,K,correlation = 'independence', Gamma=NULL,Beta,Sigma2,Rho=NULL)
{	
        N <- length(unique(id))
        m=nrow(X)
        
        Sigma <- function(d,Sigma2,Rho=NULL,chol=TRUE,correlation){
                n <- ncol(as.matrix(d))
                if(correlation == 'exponential'){ S <- Sigma2*exp(-(as.matrix(d)/Rho))}
                else if(correlation =='gaussian'){ S <- Sigma2*exp(-(as.matrix(d)/Rho)^2)}
                else if(correlation =='random slope'){ S <- (Sigma2 - Rho*Sigma2)*diag(n) + cbind(1,1:n)%*%matrix(c(Rho*Sigma2,0,0,Sigma2*.05),nrow=2,ncol=2)%*%t(cbind(1,1:n))}
                else if(correlation == 'compound'){ S <- matrix(Sigma2,nrow=n,ncol=n)*(matrix(Rho,nrow=n,ncol=n)-diag(Rho-1,n))}
                else if(correlation == 'independence'){ S <- diag(Sigma2,nrow=n,ncol=n)}
                else{ stop("Specified 'correlation' not found")}
                if(chol){
                        return(t(chol(S)))
                }else{ return(S) }
        } 
        
        z=rep(1,m)
        y=rnorm(m)
        cl = rep(0,length(y))
        for(i in unique(id)){
                l =sample(1:K,1,replace=TRUE,prob=invlogit(1,Gamma))
                y[which(id==i)]=Sigma(dist(Time[which(id==i)]),Sigma2[l],Rho[l],chol=TRUE,correlation)%*%y[which(id==i)]+as.matrix(cbind(1,X[which(id==i),]))%*%as.matrix(Beta[l,])		
                cl[which(id == i)] = l
        }
        return(data.frame(y=y,id=id,X=X,z=z,cl=cl))
}



invlogit <- function(z,Gamma){
        Mat <- c(exp(as.matrix(z)%*%t(as.matrix(Gamma))),1)
        Mat <- Mat/(sum(Mat))
        return(Mat)
}



dmvnorm = function(x,mean,sigma,log=FALSE){
        if(is.vector(x)){
                x = matrix(x, ncol = length(x))	
        }
        if(missing(mean)){
                mean = rep(0,length=ncol(x))
        }
        if(missing(sigma)){
                sigma = diag(ncol(x))
        }
        if(ncol(x) != ncol(sigma)){
                stop('x and sigma have non-conforming size')
        }
        if(!isSymmetric(sigma,tol = sqrt(.Machine$double.eps),check.attributes=FALSE)){
                stop('sigma must be a symmetric matrix')
        }	
        if(length(mean) != nrow(sigma)){
                stop('mean and sigma have non-conforming size')
        }
        logdet = log(det(sigma))
        center = x - mean
        distval = center%*%solve(sigma,t(center))
        logretval = -(ncol(x) * log(2 * pi) + logdet + distval)/2
        if(log){ return(logretval)
        }else{ return(exp(logretval))}
}

update.prob = function(Y,X,Z,Time,IDS,N,K,Gamma.iter,Beta.iter,Sigma.iter,correlation){
        pp <- matrix(0,nrow=length(IDS),ncol=K)
        S  <- vector('list',K)
        for(k in 1:K){
                if(correlation == 'expplus'){
                        S[[k]] <- lapply(IDS,function(i){ Sigma.iter[k,1]*((1-Sigma.iter[k,2])*exp(-as.matrix(dist(Time[[which(IDS == i)]]))*Sigma.iter[k,3]) +Sigma.iter[k,2])}) #incorporate covariance structure here
                }
                if(correlation == 'independence') S[[k]] <- lapply(IDS,function(i){ Sigma.iter[k,1]*diag(length(Time[[which(IDS == i)]]))}) #incorporate covariance structure here
        } 
        lik = function(i,k){
                invlogit(Z[[which(IDS == i)]],Gamma.iter)[k]*dmvnorm(Y[[which(IDS == i)]],drop(X[[which(IDS == i)]]%*%Beta.iter[k,]),S[[k]][[which(IDS == i)]])
        }
        
        for(k in 1:K){
                pp[,k] <- unlist(sapply(IDS,lik,k = k))
        }

        pp <- t(apply(pp,1,function(a) a/sum(a)))
        return(pp)
}

update.parameters = function(Y,X,Z,Time,IDS,N,K,pp,Gamma.iter,Beta.iter,Sigma.iter,correlation){ 
 
        log.lik.sigma = function(p,Y,X,Time,IDS,pp,correlation,Beta.iter){
                Sigma.iter = as.numeric(p)
                if(Sigma.iter[1] <= 0) return(NA)
                if(correlation == 'expplus'){
                        if(Sigma.iter[3] <= 0 | Sigma.iter[2]<= -.04/Sigma.iter[3] | Sigma.iter[2] >= 1) return(NA)#double check this works
                        S <- lapply(IDS,function(i){ Sigma.iter[1]*((1-Sigma.iter[2])*exp(-as.matrix(dist(Time[[which(IDS == i)]]))*Sigma.iter[3]) +Sigma.iter[2])}) 
                        #if(any(eigen(S[[1]])$values < 0) return(NA)
                }
                if(correlation == 'independence') S <- lapply(IDS,function(i){ Sigma.iter[1]*diag(length(Time[[which(IDS == i)]]))}) 
                
                if(any(sapply(S,function(v) any(eigen(v)$values < 0)))){ return(NA)
                }else{return(-sum(pp*sapply(IDS, function(i) dmvnorm(Y[[which(IDS == i)]],drop(X[[which(IDS == i)]]%*%Beta.iter),S[[which(IDS == i)]],log=TRUE))))}
        } 
        Yt = unlist(Y)
        Xt = do.call(rbind.fill.matrix,X)[,-1]
        ppt = do.call(rbind.fill.matrix,lapply(1:length(IDS),function(i) matrix(rep(pp[i,],length(Y[[i]])),ncol=ncol(pp),byrow=TRUE)))
        for(k in 1:K){
                #figure out a better optimizer with Beta.....weighted LS and then use optim?
                Beta.iter[k,] = lm(Yt~Xt,weights=ppt[,k])$coef
                if(correlation == 'independence'){
                        foo = optimize(log.lik.sigma,interval=c(0,300),Y=Y,X=X,Time=Time,IDS=IDS,pp=pp[,k],correlation=correlation,Beta.iter=Beta.iter[k,])$minimum
                }else{
                        foo = optim(Sigma.iter[k,],fn = log.lik.sigma,Y=Y,X=X,Time=Time,IDS=IDS,pp=pp[,k],correlation=correlation,Beta.iter=Beta.iter[k,])$par
                }
                Sigma.iter[k,] = foo
        }
        Ztmp <- do.call(rbind.fill.matrix,Z) #unlist the original data
        Gamma.iter <- summary(multinom(pp~Ztmp-1,trace=FALSE))$coeff
        Gamma.iter <- matrix(t(t(rbind(0,Gamma.iter))-Gamma.iter[K-1,])[-K,],nrow=K-1)
        return(list(Gamma.iter,Beta.iter,Sigma.iter))
}

start.em <- function(dat,IDS,K){ #improve this...to get a better starting value
        pp = t(sapply(IDS,function(i)  diag(K)[,sample(1:K,1)]))
        return(pp)
}


MixLik <- function(Y,X,Z,Time,IDS,dat,id,K,Gamma.iter,Beta.iter,Sigma.iter,correlation){
        S  <- vector('list',K)
        for(k in 1:K){
                if(correlation == 'expplus'){
                        S[[k]] <- lapply(IDS,function(i){ Sigma.iter[k,1]*((1-Sigma.iter[k,2])*exp(-as.matrix(dist(Time[[which(IDS == i)]]))*Sigma.iter[k,3]) + Sigma.iter[k,2])}) 
                }
                if(correlation == 'independence') S[[k]] <- lapply(IDS,function(i){ Sigma.iter[k,1]*diag(length(Time[[which(IDS == i)]]))}) 
        }        
        log.lik = function(i,K){
                log(sum(sapply(1:K,function(k) invlogit(Z[[i]],Gamma.iter)[k]*dmvnorm(Y[[which(IDS == i)]],drop(X[[which(IDS == i)]]%*%Beta.iter[k,]),S[[k]][[which(IDS == i)]]))))
        }
        
        sm = sum(sapply(1:length(IDS),log.lik,K))	
        return(sm)
}



#Correlation Structures: Generalized Exponential Plus Intercept, Independence 
#To Do: Speed up update.parameters function (limit parameter space), Improve Start_EM
EM <- function(form,dat,id,K,corForm = NULL,concomit = ~1, correlation='independence', tol=1E-15, max.iter = 50,verbose=TRUE){
        
        MOD <- model.matrix(form, dat) #y ~ x + x2
        p <- ncol(MOD)
        if(!is.null(concomit)){
                MOD2 <- model.matrix(concomit,dat)
                q <- ncol(MOD2)
        }
        IDS <- unique(id)
        Y <- lapply(IDS,function(i) model.extract(model.frame(form,dat),'response')[which(id==i)]) #least efficient part.
        X <- lapply(IDS,function(i) model.matrix(form,dat)[which(id==i),])
        Z <- lapply(IDS,function(i) matrix(model.matrix(concomit,dat)[which(id==i)[1],],ncol=q))
        Time <- lapply(IDS,function(i) model.matrix(corForm,dat)[which(id==i),-1])
        N <- length(IDS)
        
        if(correlation == 'independence') l <- 1
        if(correlation == 'expplus') l <- 3
        
        #Save space for parameters
        Beta.iter=matrix(rep(lm(form,dat)$coef,K),K,p,byrow=TRUE)
        Sigma.iter=matrix(0.5,K,l)
        if(correlation == 'expplus') Sigma.iter[,2] = 0
        Gamma.iter=matrix(0,nrow=K-1,ncol=q)
        
        #initialize each ID to random group or initialize parameters
        pp <- start.em(dat,IDS,K)
        
        #update parameters
        U <- update.parameters(Y,X,Z,Time,IDS,N,K,pp,Gamma.iter,Beta.iter,Sigma.iter,correlation)
        Gamma.iter <- U[[1]]
        Beta.iter <- U[[2]]
        Sigma.iter <- U[[3]]

        
        #update Log likelihood
        Loglik <- MixLik(Y,X,Z,Time,IDS,dat,id,K,Gamma.iter,Beta.iter,Sigma.iter,correlation)
        
        PrevLoglik <- log(0)
        iter <- 0	
        
        while(Loglik-PrevLoglik>tol & iter < max.iter){
                iter <- iter+1	
                
                #update Posterior Probabilities
                pp <- update.prob(Y,X,Z,Time,IDS,N,K,Gamma.iter,Beta.iter,Sigma.iter,correlation)
                if(any(apply(pp,2,sum)<0.9)) stop('Warning: K is too large')
                #update parameters
                U <- update.parameters(Y,X,Z,Time,IDS,N,K,pp,Gamma.iter,Beta.iter,Sigma.iter,correlation)
                Gamma.iter <- U[[1]]
                Beta.iter <- U[[2]]
                Sigma.iter <- U[[3]]

                
                
                PrevLoglik <- Loglik
                Loglik <- MixLik(Y,X,Z,Time,IDS,dat,id,K,Gamma.iter,Beta.iter,Sigma.iter,correlation)
                if(verbose==TRUE) cat('iter =',iter,'\nLL =',Loglik,'\ndiff =',Loglik-PrevLoglik,'\n\n')
        }
        j = length(Gamma.iter)+length(Beta.iter)+length(Sigma.iter)
        bic = -2*Loglik + j*log(N)
        #print(iter)
        return(list(Gamma = Gamma.iter, Beta = Beta.iter, Sigma = Sigma.iter,Loglik = Loglik,pp = pp, bic = bic))
}
        

N = 100
m = 5
K = 2
dat = gen(rep(1:m,N),as.matrix(rep(1:m,N)),rep(1:N,each=m),K=2,correlation='exponential',Gamma=0,Beta=matrix(c(1,10,3,9),nrow=2,byrow=TRUE),Sigma2=c(4,4),Rho=c(5,5))
dat$y.norm = unlist(tapply(dat$y,dat$id,function(v) v-mean(v,na.rm=TRUE)))
#pp <- matrix(0,nrow=length(unique(dat$id)),ncol=K)
#pp[,1] = unique(dat[,c('id','cl')])[,2]-1
#pp[,2] = 1-pp[,1]


#When is this going to make a difference???? High overlap
#Will it impact choice of K? yes.. fit each model 5 times and chose the one with the highest likelihood
foo1 = rep(NA,5)
foo2 = rep(NA,5)
for(k in 2:6){
        m1.bic = Inf
        m2.bic = Inf
        for(j in 1:5){
                tmp1 = EM(form = y.norm~X,dat=dat,id=dat$id,K = k,corForm = ~X,concomit = ~1, correlation = 'expplus',tol=1E-05, max.iter = 50,verbose=TRUE)
                tmp2 = EM(form = y.norm~X,dat=dat,id=dat$id,K = k,corForm = ~X,concomit = ~1, correlation = 'independence',tol=1E-05, max.iter = 50,verbose=TRUE)
                if(tmp1$bic < m1.bic) m1 = tmp1
                if(tmp2$bic < m2.bic) m2 = tmp2
        }
foo1[k-1]=m1$bic
foo2[k-1]=m2$bic
}

#Truth was 2
#ExpPlus recommended 3
#Independence recommended 4


#m3 = EM(form = y~X,dat=dat,id=dat$id,K = 3,corForm = ~X,concomit = ~1, correlation = 'expplus',tol=1E-5, max.iter = 50,verbose=TRUE)
#m4 = EM(form = y~X,dat=dat,id=dat$id,K = 3,corForm = ~X,concomit = ~1, correlation = 'independence',tol=1E-5, max.iter = 50,verbose=TRUE)

