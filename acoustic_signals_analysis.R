# Load the dataset of timestamps for vocalizations
# Each time step in the bird data represents 0.5 seconds
load("./vocalization_data.RData")

# Define the limit as to how many time steps can be saved in the memory kernel
# Theoretically, this could be infinite, but that makes model fitting difficult
# So we set it to a reasonable upper bound of 60 timesteps
max.memory <- 60

# Define a cutoff in how many total songs are held in the memory kernel
# We have set it equal to maxmem
# This means if there were 60 vocalizations in 60 timesteps, an individual could remember all of them
# We could make the cutoff more strict if there were a biological reason
# For example, if we knew how memories of neighbor vocalizations were stored by birds
cutoff <- 60

# We are fitting all 15 ovenbird sites + 1 bullfrog site separately
# Each site could be run in a separate script, but we will run them all in a loop
# If you want to run only one site, set a value for s (e.g. s <- 1) and run the code inside the loop
for (s in 1:length(vocalizations)) {
  
  # Subset the list of data to the current site
  site_data <- vocalizations[[s]]
  
  # Create an array containing the history of events leading up to each time stamp
  history <- array(dim = c(dim(site_data)[1], dim(site_data)[2], max.memory))
  
  for (i in 1:dim(site_data)[1]) {
    for (j in 1:dim(site_data)[2]) {
      prior_calls <- which(site_data[i, c(1:j)] > 0) # detect vocalizations
      if (any(prior_calls) > 0) {
        events <- j - prior_calls # current time minus the history of the process
        memories <- events[which(events < cutoff)] # history of process at each time within our constraints
      } else{
        memories <- 0 # if no calls within our constraints
      }
      length(memories) <- max.memory
      memories[is.na(memories)] <- 0
      
      history[i, j,] <- memories
    }
  }
  
# Create an object containing our data to pass to JAGS
  # note that we pass site_data twice; this is to separately estimate the models for comparison
  # 1) Self-exciting point process model
  # 2) A modified Poisson model with a variable lambda to see if the increased specificity in our model changed fit
  # 3) A standard Poisson model without a variable lambda
  jags.data <- list(
    hwk.calls = site_data,
    m.pois.calls = site_data,
    pois.calls = site_data,
    history = history,
    max.memory = max.memory,
    n.obs = dim(site_data)[2],
    n.sites = dim(site_data)[1]
  )
  rm(history)
  # Pass all our data to JAGS
  site_model <- R2jags::jags(
    data = jags.data,
    parameters.to.save = c(
      "mu",
      "sepp.lambda",
      "gamma",
      "alpha",
      "beta",
      "sepp.sim",
      "pois.lambda",
      "pois.sim",
      "m.pois.lambda",
      "m.pois.sim"
    ),
    model.file = "./scripts/submission/acoustic_signals_jags.R",
    n.chains = 3,
    n.iter = 3000,
    n.burnin = 1000,
    n.thin = 3
  )
  rm(jags.data)
  # save the model output using the name of the site
  save(site_model, file=paste("./", names(all_sites)[s], ".rda", sep=""))
}
