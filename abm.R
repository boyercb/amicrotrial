library(ABM)
library(tidyverse)
library(scico)
library(ggokabeito)


clus = 4 # number of clusters in the two-stage trial
sims = 5 # number of simulations to run per cluster,

runsim <- function (
  # parameters
  gammaR = 0.2, # the recovery rate for R (1/30)
  gammaS = 0.2, # the recovery rate for S (1/30)
  beta = 0.4, # the transmission rate
  c = 0.02, # the fitness cost of resistance
  n = 1000, # the initial population size
  mu = 0.01, # the discharge/death rate
  lambda = 10, # the admission rate
  m = 0.75, # probability of X admission (versus S)
  S0 = 5, # the initial number of infectious agents with susceptible infection
  R0 = 5, # the initial number of infectious agents with resistant infection
  sigma = 1, # the rate of symptom acquisition after infection
  rho = 0.5, # the probability of treatment given symptoms
  pi = 0.8, # the probability of receiving drug A versus drug B
  tau = 0.7, # the increase in recovery rate of covering drug
  timesteps = 200 # number of timesteps 
){
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
  sim$addLogger(newCounter("TS1_inc", "S", "TS1"))
  sim$addLogger(newCounter("TS2_inc", "S", "TS2"))
  sim$addLogger(newCounter("UR", "UR"))
  sim$addLogger(newCounter("TR1", "TR1"))
  sim$addLogger(newCounter("TR2", "TR2"))
  sim$addLogger(newCounter("TR1_inc", "R", "TR1"))
  sim$addLogger(newCounter("TR2_inc", "R", "TR2"))
  
  # use random mixing within clusters for now
  mx = newRandomMixing()
  
  # initialize contacts
  sim$addContact(mx)
  
  # susceptible infection treatment
  sim$addTransition("S"->"TS1", sigma * rho * pi)
  sim$addTransition("S"->"TS2", sigma * rho * (1 - pi))
  sim$addTransition("S"->"US", sigma * (1 - rho))
  # sim$addTransition("S"->"US", sigma, 
  #                   changed_callback = function(time, agent) {
  #                     if (runif(1) < rho) {
  #                       if (runif(1) < pi) {
  #                         setState(agent, "TS1")
  #                       } else {
  #                         setState(agent, "TS2")
  #                       }
  #                     } else {
  #                       setState(agent, "US")
  #                     }
  #                   })

  # resistant infection treatment
  sim$addTransition("R"->"TR1", sigma * rho * pi)
  sim$addTransition("R"->"TR2", sigma * rho * (1 - pi))
  sim$addTransition("R"->"UR", sigma * (1 - rho))
  # sim$addTransition("R"->"UR", sigma, 
  #                   changed_callback = function(time, agent) {
  #                     if (runif(1) < rho) {
  #                       if (runif(1) < pi) {
  #                         setState(agent, "TR1")
  #                       } else {
  #                         setState(agent, "TR2")
  #                       }
  #                     } else {
  #                       setState(agent, "UR")
  #                     }
  #                   })
  
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
  
  x = sim$run(0:timesteps)
  return(x)
}

# initialize results list
res <- list()

# repeatedly run simulation 

# loop over clusters in trial 
for (i in 1:clus) {
  # assign pi: two-arm 50/50 split
  if (i <= floor(clus / 2)) {
    pi = 0.8
  } else {
    pi = 0.5
  }
  
  # loop over simulations per cluster
  for (j in 1:sims) {
    
    # run simulation
    x = runsim(pi = pi, timesteps = 50)
    
    x$clus = i
    x$sim = j
    
    # append to list
    res[[(i - 1) * sims + j]] = x
  }
}


# show example results
print(res[[1]])

# create one big data.frame
res <- bind_rows(res, .id = "id")

# reorganize for plotting
plot_df <-
  res |>
  mutate(
    R = R + TR1 + TR2 + UR,
    S = S + TS1 + TS2 + US,
    X = X,
    sim = as.numeric(sim)
  ) |>
  select(sim, times, clus, R, S , X) |>
  pivot_longer(-c(sim, times, clus))

plot_df2 <-
  res |>
  mutate(
    drugA = TR1_inc + TS1_inc,
    drugB = TR2_inc + TS2_inc,
    sim = as.numeric(sim)
  ) |>
  select(sim, times, clus, drugA, drugB) |>
  pivot_longer(-c(sim, times, clus))

# create some plots
ggplot(plot_df, aes(x = times, y = value, color = name)) +
  facet_grid(clus~sim, labeller = label_both) +
  geom_line(linewidth=1) +
  scale_color_scico_d(palette = "batlow") +
  theme_minimal()

ggplot(plot_df2, aes(x = times, y = value, color = name)) +
  facet_wrap(~sim, labeller = label_both) +
  geom_line(linewidth=1) +
  scale_color_okabe_ito() +
  theme_minimal()

plot_df3 <-
  res |>
  mutate(
    drugA = (TR1_inc + TS1_inc) / (TR1_inc + TS1_inc + TR2_inc + TS2_inc),
    drugB = (TR2_inc + TS2_inc) / (TR1_inc + TS1_inc + TR2_inc + TS2_inc),
    sim = as.numeric(sim)
  ) |>
  select(sim, times, clus, drugA, drugB) |>
  pivot_longer(-c(sim, times, clus))

ggplot(plot_df3, aes(x = times, y = value, color = name)) +
  facet_grid(clus~sim, labeller = label_both) +
  geom_line(linewidth=1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_color_okabe_ito() +
  theme_minimal()
