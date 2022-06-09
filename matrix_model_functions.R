## Matthew Coates


## IGNORES MIGRATION
## function needs to do male and female at the same time so that births can happen correctly
sim_markov_chain_array_sex <- function(x0,pf,pm,n_cycles,xadd,state_names,base_year,bprev,b_enroll_b,b_enroll_d){
  ## x0 is initializing pop
  ## pf and pm are transition probability matrices (array) for males and females
  ## n_cycles is number of cycles to run forward
  ## xadd is fertility to calculate births
  ## state_names is state names
  ## base_year is base year of running model
  ## bprev is birth prevalence of relevant diseases if needed (e.g. for SCD)
  ## b_enroll_b is baseline % of births with SCD enrolled immediately in care
  ## b_enroll_d is the % point change in this for each year
  
  ## initialize
  ages <- length(unique(x0$age))
  dimensions <- c(ages,length(state_names),n_cycles) ## (dimensions should be states x ages x cycles)
  male <- array(NA, dim=dimensions, dimnames=list(as.character(0:(length(unique(x0$age))-1)),state_names,1:n_cycles)) # Initialize Markov trace
  female <- copy(male)
  
  ## append on the starting populations
  male <- abind(as.matrix(x0[sex=="male",state_names,with=F]), male, along=3) # Markov trace at cycle 0 is initial state, but we have initial state for 
  female <- abind(as.matrix(x0[sex=="female",state_names,with=F]), female, along=3) # Markov trace at cycle 0 is initial state, but we have initial state for 
  
  ## initialize prevalence outputs (correction to get "mid-cycle" populations)
  pmale <- copy(male)
  pfemale <- copy(female)
  
  ## ensure that the length of the transition matrix is the same as the number of cycles*age
  if (!dim(pm)[3]==(n_cycles+1)*ages) stop("issue with dimensions of transition array")

  ## loop over cycles
  for (t in 1:n_cycles){ # Simulating state vectors at each cycle with for loop
    
    b_enroll <- b_enroll_b + b_enroll_d*t ## proportion of births with SCD enrolled in care immediately through newborn screening
    if (b_enroll > 0.95) b_enroll <- 0.95
    
    y <- base_year+t
    
    ## first, estimate the number of new births there should be based on fertility (test more efficient way of doing this with fewer steps to save time)
    births <- sum((apply(female[,1:4,t],MARGIN=1,sum)[11:55])*xadd[year==y-1]$val)
    malebirths <- 1.05/2.05*births
    femalebirths <- births-malebirths
    
    ## susceptible, prevalent, enrolled, ltfu, dead, dead
    maleentry <- c(malebirths*(1-bprev[sex=="male"]$p_prev),malebirths*bprev[sex=="male"]$p_prev*(1-b_enroll),malebirths*bprev[sex=="male"]$p_prev*(b_enroll),0,0,0)
    femaleentry <- c(femalebirths*(1-bprev[sex=="female"]$p_prev),femalebirths*bprev[sex=="female"]$p_prev*(1-b_enroll),femalebirths*bprev[sex=="female"]$p_prev*(b_enroll),0,0,0)
    
    ## age population forward a year and add births
    mstart <- rbind(maleentry,male[,,t][1:(nrow(male[,,t])-1),])
    fstart <- rbind(femaleentry,female[,,t][1:(nrow(female[,,t])-1),])
    
    ## then, apply matrix multiplication to get new pops in each state
    for (ageindex in 1:nrow(male[,,t + 1])) {
      ## we can find the slice of the transition matrix to use with this equation: (t-1)*ages+ageindex
      male[ageindex,,t + 1] <- mstart[ageindex,] %*% pm[,,(t-1)*ages+ageindex]
      female[ageindex,,t + 1] <- fstart[ageindex,] %*% pf[,,(t-1)*ages+ageindex]
      
    }
    ## calculate prevalence as midpoint of cycle (for output)
    pmale[,,t + 1] <- (male[,,t + 1]+mstart)/2
    pfemale[,,t + 1] <- (female[,,t + 1]+fstart)/2
  }
  return(list(pmale,pfemale))
}












