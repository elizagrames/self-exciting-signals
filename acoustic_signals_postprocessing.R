# This script is for processing the output of the self-exciting point process
# models and Poisson models for vocalization data

# Functions --------------------------------------------------------------------
# Get credible intervals on differences between model estimates and observed values
ci_counts <- function(mod, obs.data){
  nsims <- dim(mod)[1]
  nsites <- dim(mod)[2]
  
  tmp.mod <- array(dim=c(nsims, nsites))
  for(s in 1:nsims){
    for(i in 1:nsites){
      tmp.mod[s,i] <- length(which(mod[s,i,]>0)) - sum(obs.data[i,])
    }
  }
  
  mod.mean <- apply(tmp.mod, 2, mean)
  mod.ci <- apply(tmp.mod, 2, quantile, c(0.025, 0.975))
  counts <- rowSums(obs.data)
  tbl <- cbind(counts, mod.ci[1,], mod.mean, mod.ci[2,])
  colnames(tbl) <- c("observed", "lower_ci", "mean", "upper_ci")
  output <- list(tbl, tmp.mod)
  return(output)
}

# Get credible intervals on the proportion of overall vocalization rate
# attributable to the self-excitement process
ci_condit <- function(param, lambda){
  nsims <- dim(param)[1]
  nsites <- dim(param)[2]
  
  tmp <- array(dim=c(nsims, nsites))

    for(s in 1:nsims){
      for(i in 1:nsites){
        if(mean(param[s,i,])!=0){
          tmp[s,i] <- mean(param[s,i,])/mean(lambda[s,i,])
        } else{tmp[s,i] <- NA}
      }
    }
    
  estimates <- apply(tmp, 2, mean)
  gest.ci <- apply(tmp, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
  tbl <- as.data.frame(cbind(gest.ci[1,], estimates, gest.ci[2,]))
  colnames(tbl) <- c("lower_ci", "mean", "upper_ci")
  output <- list(tbl, tmp)
  return(output)
}

# Function to nicely plot posterior densities
plot_density <-
  function(input, xmin, xmax, color, ymin, ymax, xlab) {
    lightcol <- paste(color, "60", sep = "")
    darkcol <- color
    plot(
      density(input[, 1]),
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      main = "",
      type = "n",
      xlab = xlab,
      axes = FALSE
    )
    axis(2)
    for (i in 1:dim(input)[2]) {
      dens <- density(input[, i])
      polygon(
        x = c(rev(dens$x), dens$x),
        y = c(rep(0, length(dens$x)), dens$y),
        col = lightcol,
        border = lightcol
      )
    }
    
    lines(density(input),
          lwd = 3,
          lty = 1,
          col = darkcol)
    abline(
      v = 0,
      lwd = 3,
      lty = 2,
      col = "black"
    )
    legend(
      "topright",
      legend = c("Point-level  ",
                 "No difference  ",
                 "Overall"),
      bty = "n",
      cex = .75,
      col = c(lightcol, "black", darkcol),
      pch = c(15, NA, NA),
      lty = c(NA, 2, 1),
      lwd = 2
    )
  }

# Generate credible intervals
get_ci <- function(model, nsims, x, y, length.out, bandwidth) {
  newx <- seq(min(x), max(x), length.out = length.out)
  pred <- array(dim = c(length(newx), nsims))
  for (i in 1:length(newx)) {
    for (j in 1:nsims) {
      pred[i, j] <- sample(model$BUGSoutput$sims.list$a, 1) + 
        sample(model$BUGSoutput$sims.list$b, 1) * newx[i]
    }
  }
  newpred <- apply(model$BUGSoutput$sims.list$predictor, 2, mean)
  predline <- cbind(exp(newpred), x)
  predline <- predline[order(predline[, 2]), ]
  
  predci <- apply(exp(pred), 1, quantile, c(0.05, 0.95))
  li <- ksmooth(newx, predci[1, ], kernel = "normal", bandwidth = bandwidth)
  ui <- ksmooth(newx, predci[2, ], kernel = "normal", bandwidth = bandwidth)
  
  output <- list(newx, predline, li$y, ui$y)
  return(output)
}

# Function to fit parameter estimates to site metadata
fit_glms <- function() {
  for (i in 1:n.obs) {
    y[i] ~ dgamma(shape, shape / exp(predictor[i]))
    predictor[i] <- a + b * x[i]
  }
  # Priors
  a ~ dnorm(0, 0.0001)
  b ~ dnorm(0, 0.0001)
  shape ~ dunif(0, 500)
  
}

# Read in files and get model estimates ----------------------------------------

# Rather than uploading all the raw output for all the models (it's a lot!) 
# that you get from running the models in JAGS, here are the summaries
load("./model_summaries.rda")
load("./vocalization_data.rda")

metadat <- read.csv("./sample_metadata.csv")

# Parameters in response to covariates -----------------------------------------

# get site-specific parameter estimates
gammas <- mus <- betas <- alphas <- samplesize <- c() # placeholders

for(i in 1:length(model_summaries)){
  tmp <- model_summaries[[i]]
  alphas[i] <- tmp[tmp$X == "alpha", 'mean']
  betas[i] <- tmp[tmp$X == "beta", 'mean']
  mus[i] <- tmp[tmp$X == "mu", 'mean']
  gammas[i] <-
    mean(tmp[tm::removePunctuation(tm::removeNumbers(tmp$X)) == "gamma", 'mean'])
  samplesize[i] <-
    length(which(tm::removeNumbers(tm::removePunctuation(tmp$X)) == "mpoislambda"))
}  



size.class <- metadat$Size..ha.
size.class[size.class<200] <- 1
size.class[size.class>200] <- 2

dates <- as.numeric(metadat$Date)

mod_mu_sm <- R2jags::jags(
  data = list(
    n.obs = length(mus[!is.na(dates)  & size.class!=2 ]),
    y = mus[!is.na(dates) & size.class!=2],
    x = dates[!is.na(dates)  & size.class!=2]
  ),
  parameters.to.save = c("a", "b", "shape", "x", "predictor"),
  model.file = fit_glms
)


mod_mu_bg <- R2jags::jags(
  data = list(
    n.obs = length(mus[!is.na(dates)  & size.class==2]),
    y = mus[!is.na(dates) & size.class==2],
    x = dates[!is.na(dates)  & size.class==2]
  ),
  parameters.to.save = c("a", "b", "shape", "x", "predictor"),
  model.file = fit_glms
)

ci_mu_bg <-
  get_ci(
    mod_mu_bg,
    nsims = 5000,
    x = dates[!is.na(dates) & size.class==2],
    y = mus[!is.na(dates) & size.class==2],
    length.out = 1000,
    bandwidth = 0.25
  )


ci_mu_sm <-
  get_ci(
    mod_mu_sm,
    nsims = 5000,
    x = dates[!is.na(dates) & size.class!=2],
    y = mus[!is.na(dates) & size.class!=2],
    length.out = 1000,
    bandwidth = 0.25
  )


plot(mus ~ dates, cex=2.5, ylim=c(0,0.08), axes=F, ylab="Background rate (mu)", 
     xlab="Date", type="n")

polygon(
  x = c(ci_mu_bg[[1]], rev(ci_mu_bg[[1]])),
  y = c(ci_mu_bg[[3]], rev(ci_mu_bg[[4]])),
  col = "#1e441e80",
  border = "#1e441e00"
)

polygon(
  x = c(ci_mu_sm[[1]], rev(ci_mu_sm[[1]])),
  y = c(ci_mu_sm[[3]], rev(ci_mu_sm[[4]])),
  col = "#4c9f7080",
  border = "#4c9f7000"
)
lines(ci_mu_sm[[2]][,1] ~ ci_mu_sm[[2]][,2], lwd=2, lty=2, col="#4c9f70")
lines(ci_mu_bg[[2]][,1] ~ ci_mu_bg[[2]][,2], lwd=2, lty=2, col="#1e441e")

points(mus ~ dates, cex=3, ylim=c(0,0.08), pch=21,
       bg=c("#4c9f70", "#1e441e")[factor(size.class==2)], col="grey10")

legend("topright", legend=c("Forests <200 ha", "Forests >200 ha", "90% CI"), 
       bty="n", pch=c(21, 21, 15), 
       pt.bg=c("#4c9f70", "#1e441e"), col=c("grey10", "grey10", "#00000060"), pt.cex=2)
par(las=1)
axis(1, cex.axis=1.5, at=seq(18, 54, 7), 
     labels = c("May 18", "June 25", "June 01", "June 08", "June 15", "June 22"))
axis(2, cex.axis=1.5)

# Example posterior density plots for ovenbirds --------------------------------

# We only include one of the model files here because each model is >2Gb 

# Note: this file is from the model run for the paper
# Based on a reviewer suggestion, the parameter names are different from the 
# contained scripts, however, the model structure is unchanged (e.g. sim_H is
# now labeled sepp.sim and sim_M is now m.pois.sim)

load("./ovenbird_example_mod.rda")
obs.data <- vocalizations[[2]]

# differences between observed and estimated data
sepp.calls <- ci_counts(site_model$BUGSoutput$sims.list$sim_H, obs.data)
m.pois.calls <- ci_counts(site_model$BUGSoutput$sims.list$sim_M, obs.data)
pois.calls <- ci_counts(site_model$BUGSoutput$sims.list$sim_P, obs.data)
    
# proportion of rate attributable to the self-excitement process
sepp.gamma.prop <- ci_condit(param=site_model$BUGSoutput$sims.list$gamma, 
                          site_model$BUGSoutput$sims.list$lambdaH)

output <- list(sepp.calls, m.pois.calls, pois.calls, sepp.gamma.prop)
names(output) <- c("sepp.calls", "m.pois.calls", "pois.calls", "sepp.gamma.prop")

xmax <- ceiling(max(unlist(lapply(output[1:3][[2]], max))) / 5) * 5
xmin <- floor(min(unlist(lapply(output[1:3][[2]], min))) / 5) * 5

par(
  mfrow = c(2, 1),
  las = 1,
  pty = "s",
  mar = c(5, 5, 1, 1)
)

plot_density(
  output[[1]][[2]],
  xmin = xmin,
  xmax = xmax,
  color = "#1282A2",
  ymin = 0,
  ymax = 0.1,
  xlab = ""
)
legend("topleft", "a", cex = 1.5, bty = "n")
plot_density(
  output[[3]][[2]],
  xmin = xmin,
  xmax = xmax,
  color = "#1282A2",
  ymin = 0,
  ymax = 0.1,
  xlab = "Estimated - observed calls"
)
legend("topleft", "b", cex = 1.5, bty = "n")
axis(1)



