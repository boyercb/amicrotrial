library(ABM)
library(tidyverse)
library(scico)
library(ggokabeito)

runsim <- function () {
  # parameters
  gammaR = 0.2 # the recovery rate for R (1/30)
  gammaS = 0.2 # the recovery rate for S (1/30)
  beta = 0.4 # the transmission rate
  c = 0.02 # the fitness cost of resistance
  n = 10000 # the initial population size
  mu = 0.01 # the discharge/death rate
  lambda = 100 # the admission rate
  m = 0.75 # probability of X admission (versus S)
  S0 = 5 # the initial number of infectious agents with susceptible infection
  R0 = 5 # the initial number of infectious agents with resistant infection
  sigma = 1 # the rate of symptom acquisition after infection
  rho = 0.5 # the probability of treatment given symptoms
  p = 0.8 # the probability of receiving drug A versus drug B
  tau = 0.1 # the increase in recovery rate of covering drug
  
  sims = 10 # number of simulations to run
  
  # a function to create an admission event, the argument "time" is the current simulation time
  admit.event = function(time) {
    newEvent(time + rexp(1, lambda), function(time, sim, agent) {
      # assign X with probability p otherwise assign S
      if (runif(1) < m) { 
        a = newAgent("X", death_time = time + rexp(1, mu))
      } else { 
        a = newAgent("S", death_time = time + rexp(1, mu))
      }
      addAgent(sim, a)
      schedule(sim, admit.event(time))
    })
  }
  
  # initialize simulation, create the agents, and seed infections
  sim = Simulation$new(as.list(c(rep("S", S0), rep("R", R0), rep("X", n-S0-R0))))
  
  # set the discharge times
  for (i in 1:n) {
    a = sim$agent(i)
    setDeathTime(a, rexp(1, mu))
  }
  
  # schedule a birth event
  sim$schedule(admit.event(0))
  
  # set up compartment loggers
  sim$addLogger(newCounter("X", "X"))
  sim$addLogger(newCounter("S", "S"))
  sim$addLogger(newCounter("R", "R"))
  sim$addLogger(newCounter("US", "US"))
  sim$addLogger(newCounter("TS1", "TS1"))
  sim$addLogger(newCounter("TS2", "TS2"))
  sim$addLogger(newCounter("UR", "UR"))
  sim$addLogger(newCounter("TR1", "TR1"))
  sim$addLogger(newCounter("TR2", "TR2"))
  
  # use random mixing within clusters for now
  mx = newRandomMixing()
  
  # initialize contacts
  sim$addContact(mx)
  
  # susceptible infection treatment
  sim$addTransition("S"->"TS1", sigma * rho * p)
  sim$addTransition("S"->"TS2", sigma * rho * (1 - p))
  sim$addTransition("S"->"US", sigma * (1 - rho))
  
  # resistant infection treatment
  sim$addTransition("R"->"TR1", sigma * rho * p)
  sim$addTransition("R"->"TR2", sigma * rho * (1 - p))
  sim$addTransition("R"->"UR", sigma * (1 - rho))
  
  # susceptible infection recovery
  sim$addTransition("TS1"->"X", gammaS + tau)
  sim$addTransition("TS2"->"X", gammaS + tau)
  sim$addTransition("US"->"X", gammaS)
  
  # resistant infection recovery
  sim$addTransition("TR1"->"X", gammaR)
  sim$addTransition("TR2"->"X", gammaR + tau)
  sim$addTransition("UR"->"X", gammaR)
  
  # new susceptible infections
  sim$addTransition("S" + "X" -> "S" + "S" ~ mx, beta)
  sim$addTransition("TS1" + "X" -> "TS1" + "S" ~ mx, beta)
  sim$addTransition("TS2" + "X" -> "TS2" + "S" ~ mx, beta)
  sim$addTransition("US" + "X" -> "US" + "S" ~ mx, beta)
  
  # new resistant infections
  sim$addTransition("R" + "X" -> "R" + "R" ~ mx, beta * (1 - c))
  sim$addTransition("TR1" + "X" -> "TR1" + "R" ~ mx, beta * (1 - c))
  sim$addTransition("TR2" + "X" -> "TR2" + "R" ~ mx, beta * (1 - c))
  sim$addTransition("UR" + "X" -> "UR" + "R" ~ mx, beta * (1 - c))
  
  x = sim$run(0:200)
  return(x)
}

# initialize results list
res <- list()

# repeatedly run simluation 
for (i in 1:sims) {
  # run simulation
  x = runsim()
  
  # append to list
  res[[i]] = x
}

# show example results
print(res[[1]])

# create one big data.frame
res <- bind_rows(res, .id = "sim")

# reorganize for plotting
plot_df <-
  res |>
  mutate(
    R = R + TR1 + TR2 + UR,
    S = S + TS1 + TS2 + US,
    X = X,
    sim = as.numeric(sim)
  ) |>
  select(sim, times, R, S , X) |>
  pivot_longer(-c(sim, times))

plot_df2 <-
  res |>
  mutate(
    drugA = TR1 + TS1,
    drugB = TR2 + TS2,
    sim = as.numeric(sim)
  ) |>
  select(sim, times, drugA, drugB) |>
  pivot_longer(-c(sim, times))

# create some plots
ggplot(plot_df, aes(x = times, y = value, color = name)) +
  facet_wrap(~sim, labeller = label_both) +
  geom_line(linewidth=1) +
  scale_color_scico_d(palette = "batlow") +
  theme_minimal()

ggplot(plot_df2, aes(x = times, y = value, color = name)) +
  facet_wrap(~sim, labeller = label_both) +
  geom_line(linewidth=1) +
  scale_color_okabe_ito() +
  theme_minimal()

