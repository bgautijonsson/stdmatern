library(stdmatern)
library(posterior)
library(tidyverse)

dim <- 20
rho_true <- 0.5
nu <- 0
n_obs <- 10

y <- rmatern_copula(n_obs, dim, rho_true, nu)

log_likelihood <- function(par, data, ...) {

  sum(dmatern_copula(data, dim, plogis(par), nu))
}

log_prior <- function(par) dnorm(par)

fn <- function(rho) {
  sum(dmatern_copula_eigen(y, dim, dim, rho, rho, nu))
}

optimize(
  log_likelihood,
  lower = 0,
  upper = 1,
  maximum = TRUE,
  data = y
)

metropolis_hastings_rho <- function(n_iterations, initial_rho, data, log_likelihood, log_prior, proposal_sd) {
  current_rho <- initial_rho
  samples <- numeric(n_iterations)
  
  for (i in 1:n_iterations) {
    # Propose a new value
    proposed_rho <- rnorm(1, mean = current_rho, sd = proposal_sd)
    
    # Calculate acceptance ratio
    log_r <- (log_likelihood(proposed_rho, data) + log_prior(proposed_rho)) - 
             (log_likelihood(current_rho, data) + log_prior(current_rho))
    
    # Accept or reject
    if (log(runif(1)) < log_r) {
      current_rho <- proposed_rho
    }
    
    # Store the sample
    samples[i] <- current_rho
  }
  
  return(samples)
}

n_samps <- 1e4
results <- metropolis_hastings_rho(
  n_samps,
  -10,
  data = y,
  log_likelihood = log_likelihood,
  log_prior = log_prior,
  proposal_sd = 0.3
)

tibble(
  logit_rho = results,
  rho = plogis(results)
) |> 
  filter(row_number() > (n() / 2)) |> 
  as_draws_df() |> 
  summarise_draws()

plot(results, type = "l")
