require(deSolve)

sir_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {  # tell R to look for variable names within the state and parameters objects
    
    N <- S+I+R
    
    # New: defining lambda as a function of beta and I:
    lambda <- beta * I/N
    
    # The differential equations
    dS <- -lambda * S               # people move out of (-) the S compartment at a rate lambda (force of infection)
    dI <- lambda * S - gamma * I    # people move into (+) the I compartment from S at a rate lambda, 
    # and move out of (-) the I compartment at a rate gamma (recovery)
    dR <- gamma * I                 # people move into (+) the R compartment from I at a rate gamma
    
    # Return the number of people in the S, I and R compartments at each timestep 
    # (in the same order as the input state variables)
    return(list(c(dS, dI, dR))) 
  })
  
}

run_sir_model <- function(beta, gamma, S0 = 999999, I0 = 1, R0 = 1, duration = 60) {
  
  initial_state_values <- c(S = S0,  # nearly the whole population we are modelling is susceptible to infection
                            I = I0,       # the epidemic starts with a single infected person
                            R = R0)       # there is no prior immunity in the population
  
  parameters <- c(beta = beta, # the force of infection, which acts on susceptibles
                  gamma = gamma)   # the rate of recovery, which acts on those infected

  times <- seq(from = 0, to = duration, by = 1) 
  
  output <- as.data.frame(ode(y = initial_state_values, 
                              times = times, 
                              func = sir_model,
                              parms = parameters))
  
  # Plotting the output
  plot(x = output$time,                                 # time on the x axis
       y = output$S,                                    # the number of susceptible people at each timestep on the y axis
       type = "l",                                      # type = "l" tells R we want lines rather than points
       ylim = c(0,(S0+I0+R0)),                             # the limits of the y axis (from 0 to the total number of people)
       xlab = "Time (days)", ylab = "Number of people") # add labels
  lines(x = output$time, y = output$I, col = "red")     # add the number of infected people at each timestep on the y axis
  lines(x = output$time, y = output$R, col = "blue")    # add the number of recovered people at each timestep on the y axis
  legend(x = "right",                                   # add a legend on the right side of the plot
         legend = c("S", "I", "R"),                     # assigning labels S, I and R
         col = c("black", "red", "blue"),               # to the black, red and blue lines respectively
         lty = c(1,1))                                  # both lines are of solid linetype (lty = 1)
  title(main = paste("beta =", beta, 
                     ", gamma =", gamma))
}