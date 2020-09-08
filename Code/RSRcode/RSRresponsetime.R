# Note: This uses the Guan & Haran (2018) see https://doi.org/10.1080/10618600.2018.1425625
# As an alternative to this method, consider using sparse.sglmm in the ngspatial package

library(dplyr)
library(splines)
library(ggplot2)

# travel_times dataframe includes travel time to an incident, lat + long, apparatus type responding,
#   station ID, distance from station, and time of incident (hour, month, year); data provided by ACFD

# load C and R code for Reduced rank Restricted Spatial Regression
source(file = "spatialPoisson.R")
dyn.load("phiCfast.so")

# specify degree and #knots for temporal spline effects
degree <- 5; knots <- 5; df <- knots + degree

formula=60*travel ~ distance + fapparatus + fstation +
  pbs(hour_frac,degree=degree,df=df,Boundary.knots=c(0,24)) + pbs(month_frac,degree=degree,df=df,Boundary.knots=c(0,12)) + bs(year_frac,degree=degree,df=df)
data=travel_times
coords=cbind(travel_times$longitude,travel_times$latitude)

# fit Poisson SGLMM
linmod <- glm(formula, data = travel_times, family=poisson(link="log"),contrasts = list(fapparatus = "contr.sum", fstation = "contr.sum")) # find starting values
starting <- list("beta"=coef(linmod),"s2"=2,"phi"=0.5) # set starting values
tuning   <- list("beta"=c(sqrt(diag(vcov(linmod)))),"s2"=0.1,"phi"=0.01,"w"=0.1) # set tuning parameters
priors   <- list("beta.normal"=c(100),"s2.IG"=c(2,2),"phi.Unif"=c(0.01, 1.5)) # set priors
adapt <- c("batchlength"=1e4,"n.batch"=100) # set Monte Carlo sample size = batchlength*n.batch

data=travel_times
coords2=cbind(travel_times$longitude,travel_times$latitude)
Y2 = as.matrix(travel_times$travel*60)
X2 <- model.matrix( ~ distance + fapparatus + fstation +
                      pbs(hour_frac,degree=degree,df=df,Boundary.knots=c(0,24)) + pbs(month_frac,degree=degree,df=df,Boundary.knots=c(0,12)) + bs(year_frac,degree=degree,df=df),
                    data=travel_times, contrasts = list(fapparatus = "contr.sum", fstation = "contr.sum"))
covfn <- covfndef(2.5) # define the covariance function, nu (2.5) is subject to change here and below

output <- poi_gp_arrpfit(O=Y2,coords=coords2,X=X2,covfn=covfn,adapt=adapt,nu=2.5,
                         starting=starting,tuning=tuning,priors=priors,rank=50,core=6,mul=2)

# model effects 
beta.samples <- output$arrp.params

beta.ci <- data.frame(name=colnames(X2),
                      mean=apply(beta.samples,2,mean),
                      lower=apply(beta.samples,2,function(x){quantile(x,.025)}),
                      upper=apply(beta.samples,2,function(x){quantile(x,.975)})
)

# example: effect of Fire Station on response time
station.coef <- beta.ci[substr(beta.ci$name,1,4)=="fsta",]
station.coef$names <- 101:110
station.coef$names <- factor(station.coef$names,levels=station.coef$names)

station_plot <- ggplot(data=station.coef, aes(x=names,y=exp(mean),ymin=exp(lower),ymax=exp(upper))) +
  theme_minimal() +
  geom_errorbar(color="red") +
  geom_point(size=2) +
  geom_hline(yintercept=1) +
  xlab("Station") + ylab("Effect Size") +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
station_plot

# example: effect of hour on respone time (fit using periodic basis spline)
hour.samples <- as.matrix(beta.samples[,16:25])
hour.basis <- as.data.frame( cbind(hour=travel_times$hour_frac,pbs(travel_times$hour_frac,degree=degree,df=df)) ) %>% arrange(hour)
hour.lines <- as.matrix(hour.basis[,2:(df+1)]) %*% t(hour.samples)

hour.basis$est <- apply(hour.lines,1,function(x){quantile(x,probs=.5)})
hour.basis$lower <- apply(hour.lines,1,function(x){quantile(x,probs=.025)})
hour.basis$upper <- apply(hour.lines,1,function(x){quantile(x,probs=.975)})

plot(hour.basis$hour,hour.basis$est,type="l")
lines(hour.basis$hour,hour.basis$lower)
lines(hour.basis$hour,hour.basis$upper)


