## Matthew Coates
## functions to prep data for matrix implementation of SCD model


paperdir <- "/papers/PEN-Plus_investment_case"
library(data.table)
library(openxlsx)
library(abind)

################################
## prep transition probabilities
################################

## all-cause mortality (can make cause-deleted)
prep_susc_dead <- function(clist=NA) {
  ## based on UN WPP mortality projection trends and GBD all-cause mortality in 2019
  ## (to preserve relationship between GBD epi estimates and mortality but use UN WPP projection patterns)
  
  ## load all-cause mortality (some other pre-processed info in here that may or may not be used)
  mortrate <- readRDS(paste0("inputs/gbd_ests.RDS"))
  mortrate <- mortrate[!age_name %in% c("Early Neonatal","Late Neonatal","Post Neonatal") & sex_name %in% c("male","female")]
  mortrate <- unique(mortrate[,c("sex_name","age","year","location_name","age_id","age_name","mx_allcauses"),with=F])
  setnames(mortrate,"sex_name","sex")
  
  ## load wpp projections
  proj <- readRDS(paste0("inputs/future_mort_scalars.RDS"))
  proj <- dcast.data.table(proj,location_name+sex+age_group_start~year,value.var = "ratio")
  proj[,`2019`:=NULL]
  mortrate[,age_group_start:=floor(age/5)*5]
  mortrate[age %in% c(1:4),age_group_start:=1]
  
  ## create projected mx_all cause
  ## in future, can change so that mx_cd instead
  mortrate <- merge(mortrate,proj,by=c("location_name","sex","age_group_start"),all.x=T)
  for (i in c(2020:2090)) {
    mortrate[[paste0(i)]] <- mortrate[[paste0(i)]]*mortrate$mx_allcauses
  }
  setnames(mortrate,"mx_allcauses","2019")
  mortrate[,year:=NULL]
  
  mortrate <- melt(mortrate,id.vars=c("location_name","sex","age_group_start","age","age_name","age_id"),variable.name="year",value.name="mx_ac")
  mortrate[,year:=as.numeric(as.character(year))]
  mortrate[,qx_ac:=1-exp(-mx_ac)]
  mortrate <- mortrate[location_name %in% clist$location_name]

  return(mortrate)
  
}


## prep with-cause mortality, treated and untreated (prep in same function since they are related)
prep_prev_dead <- function(acmort=NA,rr=NA,clist=NA) {
  ## requires all-cause mortality (acmort)
  ## rr to account for effect of services
  ## clist to manage countries, make sure no mergeing issues, etc.
  
  ######################################################################
  ## without treatment
  ######################################################################
  
  ## sickle:
  sickle <- readRDS(paste0("./inputs/sickle_excess.RDS"))
  sickle <- sickle[location_name %in% clist$location_name]
  sickle[,c("lower","upper"):=NULL]
  sickle[,cause_name:="Sickle cell disorders"]
  sickle <- merge(sickle, acmort,by=c("age","sex","location_name"),all=T)
  sickle[,qx_wc_u:=qx_ac+qx_excess] ## return to this (in terms of qx versus mx in calculating excess initially)
  

  ######################################################################
  ## load effect sizes to get the with-cause mortality among ppl treated
  ######################################################################
  # eff <- data.table(read.xlsx(paste0("inputs/placeholder_effect_inputs.xlsx"))) ## for other conditions once we add
  
  ## sickle:
  sickle[,qx_wc_t:=qx_ac+qx_excess*rr]
  
  dismort <- rbind(sickle) ## combine across conditions
  if (any(c(dismort$qx_wc_u,dismort$qx_wc_t) > 1 | c(dismort$qx_wc_u,dismort$qx_wc_t) < 0)) stop("Incorrect probability value")
  dismort <- dismort[,c("age","sex","location_name","year","cause_name","qx_ac","qx_wc_u","qx_wc_t"),with=F]
  return(dismort)
  
}


prep_susc_prev <- function(clist) {
  
  ## for sickle cell, prevalence is from birth, which is accounted for in prep_birth_cohorts function
  ## here, transition probability will be 0 across ages
  sickle <- data.table(expand.grid(location_name=clist$location_name,sex=c("female","male"),age=c(0:95),year=c(2019:2090),stringsAsFactors = F))
  sickle[,p_susc_prev:=0]
  sickle[,cause_name:="Sickle cell disorders"]
  
  p_susc_prev <- rbind(sickle)
  return(p_susc_prev)
  
}


## prep birth cohorts to be entering
prep_birth_cohorts <- function() {
  
  # if doing fertility method:
  ## load asfr scalars from WPP
  fertscale <- readRDS("./inputs/future_fert_scalars.RDS")

  ## load age-specific fertility rate from GBD
  fert <- fread(paste0(paperdir,"/inputs/gbd/IHME_GBD_2019_FERTILITY_1950_2019_ASFR_Y2020M10D27.CSV"))
  fert[grepl("Tanzania",location_name),location_name:="Tanzania"]
  fert <- fert[location_name %in% fertscale$location_name & year_id == 2019]
  fert[,age_start:=as.numeric(substr(age_group_name,1,2))]
  fert <- merge(fert[,c("location_name","age_start","val"),with=F],fertscale,by=c("location_name","age_start"),all.x=T)
  fert[,val:=val*rat]
  if (any(is.na(fert$val))) stop("issue")
  fert[,rat:=NULL]
  
  expand <- data.table(data.frame(age=c(10:54)))
  expand[,age_start:=floor(age/5)*5]
  fert <- merge(fert,expand,by=c("age_start"),all.x=T,allow.cartesian = T)
  fert[,age_start:=NULL]
  
  ## the above age-specific fertility numbers can be applied to population to find births
  ## most conditions will have 0 prevalence at birth and a negligible amount under age 1, so age 0 prevalence usually doesn't matter
  ## but for SCD at least, the prevalence should be given at birth
  sickle <- readRDS(paste0("./inputs/sickle_prob_birth.RDS"))
  sickle <- sickle[year==2019]
  sickle[,year:=NULL]
  
  births <- list(fert,sickle)
  names(births) <- c("fert","Sickle cell disorders")
  
  return(births)
  
}


## prep baseline state numbers
prep_baseline_numbers <- function(clist=NA,be=NA) {
  ## be baseline coverage
  
  ## strategy for this will be different by condition since the estimates we trust come from different places
  denom <- readRDS("./inputs/coverage_denoms.RDS")
  denom <- denom[!grepl("World Bank",location_name)]
  denom <- denom[location_name %in% clist$location_name]
  denom[,age_group_start:=as.numeric(gsub(" ","",substr(age_name,1,2)))]
  denom[age_id == 28,age_group_start:=0]
  denom[,age_group_end:=age_group_start+4]
  denom[age_id == 28,age_group_end:=0]
  denom[age_id == 5,age_group_end:=4]
  denom[age_group_start==95,age_group_end:=95]
  denom[,sex:=tolower(sex_name)]

  expand <- data.table(data.frame(age=c(0:95)))
  expand[,age_group_start:=floor(age/5)*5]
  expand[age %in% c(1:4),age_group_start:=1]
  expand[,age_group_end:=age_group_start+4]
  expand[age == 0,age_group_end:=0]
  expand[age %in% c(1:4),age_group_end:=4]
  expand[age==95,age_group_end:=95]
  
  ## get rates in single-year age groups
  denom <- merge(denom,expand,by=c("age_group_start","age_group_end"),allow.cartesian = T)
  
  ## merge on populations
  pop <- fread(paste0(paperdir,"/inputs/gbd/IHME_GBD_2019_POP_SYA_2019_Y2021M01D28.CSV"))
  pop[grepl("Tanzania",location_name),location_name:="Tanzania"]
  pop <- pop[location_id %in% unique(denom$location_id) & sex_name %in% c("male","female")]
  pop[,age:=as.numeric(age_group_name)]
  pop[age_group_id == 28,age:=0]
  pop[age_group_id == 235,age:=95]
  setnames(pop,c("val","sex_name"),c("pop","sex"))
  
  denom <- merge(denom,pop[,c("sex","location_name","age","pop"),with=F],by=c("sex","location_name","age"),all.x=T)
  denom[age > 95,pop:=0]
  if (any(is.na(denom$pop))) stop("missing pop")  
  ## replace counts by multiplying aggregate rate by single-year pop
  denom[,Number:=Rate*pop]

  ## now make baseline states for each
  ##state_names <- c("Susceptible", "PrevUnenrolled", "Enrolled", "LTFU", "DeadDis","DeadNoDis")
  denom[,Susceptible:=pop-Number]
  denom[,PrevUnenrolled:=Number*(1-be)]
  denom[,Enrolled:=Number*be]
  denom[,LTFU:=0]
  denom[,DeadDis:=0]
  denom[,DeadNoDis:=0]
  
  denom <- denom[,c("location_name","cause_name","sex","age","Susceptible", "PrevUnenrolled", "Enrolled", "LTFU", "DeadDis","DeadNoDis"),with=F]
  
  return(denom)
    
}



################################################
## Functions used in creating transition matrix/array
################################################
## function to name transition matrices
name_matrix <- function(x,state_names) {
  colnames(x) <- rownames(x) <- state_names
  return(x)
}


## make probability of death = 1 in terminal age
terminal_age <- function(p){
  ## Susceptible
  p["Susceptible",] <- c(0,0,0,0,0,1) 
  p["PrevUnenrolled",] <- c(0,0,0,0,1,0)
  p["Enrolled",] <- c(0,0,0,0,1,0)
  p["LTFU",] <- c(0,0,0,0,1,0)
  return(p)
}

## check that transition values add to 1 for each state
check_tmat <- function(x) {
  ## x is assumed to be a matrix
  if (!all(abs(apply(x,MARGIN=1,sum)-1) < .000000001)) stop("Problem with transition matrix")
}

## check that transition values add to 1 for each state (in array)
check_tarray <- function(x) {
  res <- apply(x,MARGIN=3,check_tmat)
  
  #all(abs(apply(x,MARGIN=1,sum)-1) < .000000001)
  #print(dim(x))
  #print(x)
}



prep_transition_probabilities <- function(clist=NA,state_names=NA,rr=NA,p_enroll=NA,p_enroll_d=NA,base_year=NA) {
  ## all cause mort
  p_susc_dead <- prep_susc_dead(clist=clist)
  ## use to get probabilities of dying for all states
  trans_probs <- prep_prev_dead(acmort=p_susc_dead,rr=rr,clist=clist)
  setnames(trans_probs,c("qx_ac","qx_wc_u","qx_wc_t"),c("p_susc_deadnodis","p_prev_dead","p_enroll_dead"))
  trans_probs[,p_ltfu_dead:=p_prev_dead]
  ## now susceptible to prevalent
  p_susc_prev <- prep_susc_prev(clist=clist)
  trans_probs <- merge(trans_probs,p_susc_prev,by=c("location_name","cause_name","sex","age","year"))
  ## probability of enrolling (PLACEHOLDER WHILE ASSEMBLING)
  trans_probs[,p_prev_enroll:=p_enroll+p_enroll_d*(year-base_year+1)]
  ## assuming no LTFU for now
  trans_probs[,p_enroll_ltfu:=0]
  trans_probs[,p_ltfu_enroll:=0]
  ## for now, cut off at 2050
  trans_probs <- trans_probs[year < 2051]
  
  trans_probs <- trans_probs[order(location_name,cause_name,sex,year,age)]
  
  ## set up all transition probabilities here
  trans_probs[,p_susc_susc:=1-p_susc_prev-p_susc_deadnodis]
  trans_probs[,p_susc_enroll:=0]
  trans_probs[,p_susc_ltfu:=0]
  trans_probs[,p_susc_dead:=0]
  
  trans_probs[,p_prev_prev:=1-p_prev_enroll-p_prev_dead]
  trans_probs[,p_prev_susc:=0]
  trans_probs[,p_prev_ltfu:=0]
  trans_probs[,p_prev_deadnodis:=0]
  
  trans_probs[,p_enroll_susc:=0]
  trans_probs[,p_enroll_prev:=0]
  trans_probs[,p_enroll_enroll:=1-p_enroll_ltfu-p_enroll_dead]
  trans_probs[,p_enroll_deadnodis:=0]
  
  trans_probs[,p_ltfu_susc:=0]
  trans_probs[,p_ltfu_prev:=0]
  trans_probs[,p_ltfu_ltfu:=1-p_ltfu_enroll-p_ltfu_dead]
  trans_probs[,p_ltfu_deadnodis:=0]
  
  trans_probs[,p_dead_susc:=0]
  trans_probs[,p_dead_prev:=0]
  trans_probs[,p_dead_enroll:=0]
  trans_probs[,p_dead_ltfu:=0]
  trans_probs[,p_dead_deadnodis:=0]
  trans_probs[,p_dead_dead:=1]
  
  trans_probs[,p_deadnodis_susc:=0]
  trans_probs[,p_deadnodis_prev:=0]
  trans_probs[,p_deadnodis_enroll:=0]
  trans_probs[,p_deadnodis_ltfu:=0]
  trans_probs[,p_deadnodis_deadnodis:=1]
  trans_probs[,p_deadnodis_dead:=0]
  
  ## make sure everyone dies in terminal age group
  trans_probs[age==max(trans_probs$age),c("p_susc_deadnodis"):=1]
  trans_probs[age==max(trans_probs$age),c("p_susc_susc","p_susc_prev","p_susc_enroll","p_susc_ltfu","p_susc_dead"):=0]
  trans_probs[age==max(trans_probs$age),c("p_prev_dead"):=1]
  trans_probs[age==max(trans_probs$age),c("p_prev_susc","p_prev_prev","p_prev_enroll","p_prev_ltfu","p_prev_deadnodis"):=0]
  trans_probs[age==max(trans_probs$age),c("p_enroll_dead"):=1]
  trans_probs[age==max(trans_probs$age),c("p_enroll_susc","p_enroll_prev","p_enroll_enroll","p_enroll_ltfu","p_enroll_deadnodis"):=0]
  trans_probs[age==max(trans_probs$age),c("p_ltfu_dead"):=1]
  trans_probs[age==max(trans_probs$age),c("p_ltfu_susc","p_ltfu_prev","p_ltfu_enroll","p_ltfu_ltfu","p_ltfu_deadnodis"):=0]

  
  trans_probs <- trans_probs[,c("location_name","cause_name","sex","year","age",
                                "p_susc_susc","p_susc_prev","p_susc_enroll","p_susc_ltfu","p_susc_dead","p_susc_deadnodis",
                                "p_prev_susc","p_prev_prev","p_prev_enroll","p_prev_ltfu","p_prev_dead","p_prev_deadnodis",
                                "p_enroll_susc","p_enroll_prev","p_enroll_enroll","p_enroll_ltfu","p_enroll_dead","p_enroll_deadnodis",
                                "p_ltfu_susc","p_ltfu_prev","p_ltfu_enroll","p_ltfu_ltfu","p_ltfu_dead","p_ltfu_deadnodis",
                                "p_dead_susc","p_dead_prev","p_dead_enroll","p_dead_ltfu","p_dead_dead","p_dead_deadnodis",
                                "p_deadnodis_susc","p_deadnodis_prev","p_deadnodis_enroll","p_deadnodis_ltfu","p_deadnodis_dead","p_deadnodis_deadnodis"),with=F]
  
  
  ## strategy to not have this arrangement take super long: make into 3D array but with an accompanying index to be able to subset easily for use in model
  trans_array <- array(data=unlist(trans_probs[,c("p_susc_susc","p_susc_prev","p_susc_enroll","p_susc_ltfu","p_susc_dead","p_susc_deadnodis",
                                           "p_prev_susc","p_prev_prev","p_prev_enroll","p_prev_ltfu","p_prev_dead","p_prev_deadnodis",
                                           "p_enroll_susc","p_enroll_prev","p_enroll_enroll","p_enroll_ltfu","p_enroll_dead","p_enroll_deadnodis",
                                           "p_ltfu_susc","p_ltfu_prev","p_ltfu_enroll","p_ltfu_ltfu","p_ltfu_dead","p_ltfu_deadnodis",
                                           "p_dead_susc","p_dead_prev","p_dead_enroll","p_dead_ltfu","p_dead_dead","p_dead_deadnodis",
                                           "p_deadnodis_susc","p_deadnodis_prev","p_deadnodis_enroll","p_deadnodis_ltfu","p_deadnodis_dead","p_deadnodis_deadnodis"),with=F]),
                       dim=c(nrow(trans_probs),6,6))
  
  # trans_array[1,,]
  trans_array <- aperm(trans_array,c(3,2,1))
  # trans_array[,,1]

  ## name rows and columns in transition matrices
  dimnames(trans_array)[[1]] <- state_names
  dimnames(trans_array)[[2]] <- state_names
  
  ## rather than creating this index, which actually takes a very large amount of storage, pass the data table with the indices
  #trans_probs[,index:=paste(location_name,cause_name,sex,year,age,collapse="_")]
  #dimnames(trans_array)[[3]] <- trans_probs$index
  
  ## check to make sure probabilities add to 1
  check_tarray(trans_array)
  
  return(list(trans_array,trans_probs[,c("location_name","cause_name","sex","year","age"),with=F]))
 
}


