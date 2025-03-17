### Running example for the paper Spatial autoregressive modelling of epidemiological data: geometric mean model proposal.
### The dataset has been simulate using the geometry of the Columbus shapefile from the spData library

### Load libraries
library(spData)
library(spdep)
library(INLA)
library(TruncExpFam)

### Load Columbus shapefile from the spData package (49 rows)
columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1]) 
columbus <- as(columbus, "Spatial") # convert the sf object into a Spatial object from the sp package
#plot(columbus)

### Obtain spatial weights matrices
# Contiguity of order 1
nb <- poly2nb(columbus, queen=TRUE) # Construct neighbors list 
Wc <- nb2mat(nb , style="W") # Construct neighborhood matrix (Row standardized)

# Inverse distance 
coords <- coordinates(columbus)
Eu <- as.matrix(1/dist(coords)) # Inverse of Euclidean distance between centroids
diag(Eu) <- 0
sum_rows <- rowSums(Eu)
Wid<- apply(Eu, 2, function(x) x/sum_rows) # Row standardization

### Simulate the mobility matrix as the proportion of population moving from each region to the others (daily average for the whole period)
n <- 49 # Set the matrix size
set.seed(42)
#random_values <- runif(n * n) # Generate random values between 0 and 1
random_values <- rtruncexp(n*n, a = 0, b = 0.5, rate = 46) # Generate random values between 0 and 0.5 very close to zero
matContact <- matrix(random_values, n, n) # Create the matrix
diag(matContact) <- 0

## Row standardize the matrix
sum_rows <- rowSums(matContact) 
sum_rows[sum_rows==0] <- 1
matContact.st<- apply(matContact, 2, function(x) x/sum_rows) 
M <- matContact.st

### Generate random population counts for each region
set.seed(42)  # For reproducibility
population_counts <- rpois(n, lambda=200)

### Simulate auto-correlated counts according to the mobility matrix (Gibbs sampling)
# Set true values
beta <- -0.5
rho <- 0.8
#rho <- 0.5
tau <- 5
sd <- sqrt(1/tau) 

set.seed(42)
reff <- rnorm(n=n, mean=0, sd=sd) # Generate normal random effect

m=5 # m: number of Gibbs sampling iterations

# n: number of observations on each Gibbs sampling iteration

# Initialize matrices to store the generated values 
mu <- matrix(rep(0,n*m), nrow=n)  # To store the mean
y <- matrix(rep(0,n*m), nrow=n)   # To store the counts
sum <- matrix(rep(0,n*m), nrow=n) # To store the sum of w_ij*y_j
r <- matrix(rep(0,n*m), nrow=n)   # To store the rates
log.r <- matrix(rep(0,n*m), nrow=n)   # To store the log or the rates

# First iteration is to generate from an uncorrelated Poisson
mu[,1] <- exp(log(population_counts) + beta + reff) # Obtain linear predictor (mean)
set.seed(42)
y[,1] <- rpois(n,lambda=mu[,1]) # Generate Poisson counts
#r[,1] <- y[,1]/flanders$pop # Obtain the rates
r[,1] <- (y[,1]+1)/population_counts # Obtain the rates
log.r[,1] <- log(r[,1]) # Obtain the logarithm of the rates (with y+1)

## Start Gibbs sampling loop
for (j in 2:m) {  
  # Start loop for generating n samples
  for (i in 1:n) {
    # Initialize sums
    sum1 <- 0
    sum2 <- 0
    # Sum of w_ij*y_j conditioned on previous iterations
    if(i==1){
      for (s in 1:n){sum2 = sum2 + M[i,s]*log.r[s,j-1]}
    } else {
      for (s in 1:(i-1)){sum1 = sum1 + M[i,s]*log.r[s,j]}
      for (s in i:n){sum2 = sum2 + M[i,s]*log.r[s,j-1]}
    }
    # Store generated values
    sum[i,j] <- sum1 + sum2
    mu[i,j] <- exp(log(population_counts[i]) + beta + rho*sum[i,j] + reff[i])
    y[i,j]  <- rpois(1, lambda=mu[i,j])
    r[i,j]  <- (y[i,j]+1)/population_counts[i]
    log.r[i,j] <- log(r[i,j])
  }
}


### Store the last simulated dataset in the shapefile
columbus$ID <- 1:49
columbus$pop <- population_counts
columbus$sim.y <- y[,m]
# Obtain spatial lags
columbus$lag.sim.r.m <- M%*%log.r[,m] # Mobility 
columbus$lag.sim.r.c <- Wc%*%log.r[,m] # Contiguity of order1 
columbus$lag.sim.r.id <- Wid%*%log.r[,m] # Inverse distance


### Fit the geometric spatial conditional models in INLA

# Priors for INLA
prior.fixed <- list(mean.intercept = 0, prec.intercept = 1.0E-5, mean = 0, prec = 1.0E-5)
prior.iid = c(0.0001, 0.0001)


## Fit the model using the mobility matrix
formula <- sim.y ~ 1 + offset(log(pop)) + lag.sim.r.m + f(ID, model="iid", param=prior.iid)
geomm_model_m <- inla( formula = formula, data=columbus@data, family="poisson",
                 control.fixed = prior.fixed, 
                 control.compute = list(cpo=TRUE, dic = TRUE, waic=TRUE))
s.m <- rbind(geomm_model_m$summary.fixed[,c(1,2,3,5)],geomm_model_m$summary.hyperpar[,c(1,2,3,5)])
cpo.m <- -sum(log(geomm_model_m$cpo$cpo))
ic.m <- rbind(c('DIC','WAIC','CPO'),cbind(geomm_model_m$dic$dic,geomm_model_m$waic$waic,unname(cpo.m)))

# Show results
s.m
ic.m


## Fit the model using the contiguity of order 1 spatial weights matrix
formula <- sim.y ~ 1 + offset(log(pop)) + lag.sim.r.c + f(ID, model="iid", param=prior.iid)
geomm_model_c <- inla( formula = formula, data=columbus@data, family="poisson",
                       control.fixed = prior.fixed, 
                       control.compute = list(cpo=TRUE, dic = TRUE, waic=TRUE))
s.c <- rbind(geomm_model_c$summary.fixed[,c(1,2,3,5)],geomm_model_c$summary.hyperpar[,c(1,2,3,5)])
cpo.c <- -sum(log(geomm_model_c$cpo$cpo))
ic.c <- rbind(c('DIC','WAIC','CPO'),cbind(geomm_model_c$dic$dic,geomm_model_c$waic$waic,unname(cpo.c)))

# Show results
s.c
ic.c

## Fit the model using the inverse distance spatial weights matrix
formula <- sim.y ~ 1 + offset(log(pop)) + lag.sim.r.id + f(ID, model="iid", param=prior.iid)
geomm_model_id <- inla( formula = formula, data=columbus@data, family="poisson",
                       control.fixed = prior.fixed, 
                       control.compute = list(cpo=TRUE, dic = TRUE, waic=TRUE))
s.id <- rbind(geomm_model_id$summary.fixed[,c(1,2,3,5)],geomm_model_id$summary.hyperpar[,c(1,2,3,5)])
cpo.id <- -sum(log(geomm_model_id$cpo$cpo))
ic.id <- rbind(c('DIC','WAIC','CPO'),cbind(geomm_model_id$dic$dic,geomm_model_id$waic$waic,unname(cpo.id)))

# Show results
s.id
ic.id
