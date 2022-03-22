model{
    
#### Self-exciting point process model ####
  
  # In this case, n.sites refers to how many samples/recordings were taken at a location
  # So a single "site" is actually one sampled point within a forest or at the pond
  for (i in 1:n.sites) {
    for (j in 1:n.obs) {
      # A vocalization at site i in time j comes from a Poisson process with a rate lambda i,j
      # Each time stamp has its own lambda
      sepp.calls[i, j] ~ dpois(sepp.lambda[i, j])
      
      # lambda comes from a background rate, mu, and a conditional intensity, gamma
      sepp.lambda[i, j] = mu + gamma[i, j]
      
      # gamma is based on the self-excitement parameter (alpha) and an exponential memory kernel (beta)
      # This is wrapped into an ifelse statement because when the history of the process is 0, conditional intensity must be 0
      # An individual cannot respond to a neighbor's vocalization that did not happen or that it cannot remember
      for (n in 1:max.memory) {
          gamma.tmp[i, j, n] = ifelse(history[i, j, n] == 0, 0, alpha * (exp(1) ^ (-beta * history[i, j, n])))
        }
      
      # Sum the gamma for each timestep included in the memory of the process
      gamma[i, j] = sum(gamma.tmp[i, j,])
      
    }
  }
  
  # priors
  alpha ~ dgamma(0.0001, 0.0001)T(0.0001, 1)
  beta ~ dgamma(0.0001, 0.0001)T(0.0001, 5)
  mu ~ dgamma(mu.shape, mu.rate)T(0.0001, 2) # for consistency with poisson
  mu.shape <- pow(x.mu,2)/pow(sd.mu,2)
  mu.rate <- x.mu/pow(sd.mu,2)
  x.mu ~ dunif(0,100)
  sd.mu ~ dunif(0,100)

  # Simulation
  
    for(i in 1:n.sites){
      for(j in 1:n.obs){
        sepp.sim[i,j] ~ dpois(sepp.lambda[i,j])
      }
    }
    
  
##### Modified Poisson process model #### 
  
    for (i in 1:n.sites) {
      for (j in 1:n.obs) {
        # lambda is allowed to vary by site, but not time 
        # This is to check if our model fit is due to more flexibility
        m.pois.calls[i, j] ~ dpois(m.pois.lambda[i])
      }
      m.pois.lambda[i] ~ dgamma(m.pois.lambda.shape, m.pois.lambda.rate)T(0.0001, 2)
    
    }
  
  # Priors for modified Poisson
  
  m.pois.lambda.shape <- pow(x.m.pois.lambda,2)/pow(sd.m.pois.lambda,2)
  m.pois.lambda.rate <- x.m.pois.lambda/pow(sd.m.pois.lambda,2)
  x.m.pois.lambda ~ dunif(0,100)
  sd.m.pois.lambda ~ dunif(0,100)
  
  # Simulation of modified Poisson
  
  for(i in 1:n.sites){
    for(j in 1:n.obs){
      m.pois.sim[i,j] ~ dpois(m.pois.lambda[i])
  }
  }
  
#### Regular Poisson model #### 
  
  for (i in 1:n.sites) {
    for (j in 1:n.obs) {
      pois.calls[i, j] ~ dpois(pois.lambda)
    }
    
  }
  
  # Priors for Poisson model
  # Note: the hyper-priors for pois.lambda are to be consistent with comparison models
  # They have no real meaning and are simply to use consistent priors

  pois.lambda ~ dgamma(pois.lambda.shape, pois.lambda.rate)T(0.0001, 2)
  pois.lambda.shape <- pow(x.m.pois.lambda,2)/pow(sd.m.pois.lambda,2)
  pois.lambda.rate <- x.m.pois.lambda/pow(sd.m.pois.lambda,2)
  x.pois.lambda ~ dunif(0,100)
  sd.pois.lambda ~ dunif(0,100)
  
  # Simulation of regular Poisson
  
  for(i in 1:n.sites){
    for(j in 1:n.obs){
      pois.sim[i,j] ~ dpois(pois.lambda)
    }
  }
  } 
