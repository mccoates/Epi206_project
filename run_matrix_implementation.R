## Matthew Coates
## Code to run matrix implementation of SCD model

rm(list=ls())
library(data.table)
library(openxlsx)
library(ggplot2)
library(stringr)

paperdir <- "/PEN-Plus_investment_case"
tmpdir <- "/tmp/pen_plus"
source("./model/matrix_implementation_prep_functions.R")
source("./model/matrix_model_functions.R")

## some model settings
base_year <- 2022
proj_yrs <- 20

##################################
## sensitivity analysis parameters
## running 6 scenarios
## first is baseline
## 2-4 are primary scenarios 1-3 (effect from baseline)
## 5-6 are sensitivity analysis of effect parameter
##################################

## relative risk for effect of intervention on with-cause mortality
rrs <- c(.19,.19,.19,.19,.1,.5) ## 0.19 default

## baseline enrollment proportion
be <- 0.05

## scenarios for coverage scale-up

## birth enrollment (baseline and delta per year)
b_enroll_b <- 0.01
b_enroll_d <- c(0,0.01,0.03,0.05,.01,.01)

## transition probability of enrollment (baseline and delta per year)
## for other age groups
p_enroll <- 0.02
p_enroll_d <- c(0,0.005,0.01,0.02,0.005,0.005)


####################################################
## List of locations in NCDI Poverty Network
####################################################
clist <- data.table(read.xlsx(paste0("./inputs/ncdi_pov_network_countries.xlsx")))
clist <- clist[!location_name %in% c("Chhattisgarh")]
clist[location_name=="Bihar",location_name:="India"] ## for now, using india

## some locations causing issues: just focus on countries not causing issues for now
clist <- clist[!location_name %in% c("Haiti","Bangladesh","Lesotho","Togo","South Sudan","Burundi","Mali","Pakistan","India")]

####################################################
## Model states
####################################################
state_names <- c("Susceptible", "PrevUnenrolled", "Enrolled", "LTFU", "DeadDis","DeadNoDis")

##############################################
## LOOP OVER SCENARIOS TO RUN MODEL
##############################################

res <- list()
for (scen in c(1:6)) {

  ####################################################
  ## Prep baseline prevalence
  ####################################################
  init <- prep_baseline_numbers(clist=clist,be=be)
  
  ####################################################
  ## Prep transition probabilities
  ####################################################
  trans <- prep_transition_probabilities(clist=clist,state_names=state_names,rr=rrs[scen],p_enroll=p_enroll,p_enroll_d=p_enroll_d[scen],base_year)
  trans_index <- trans[[2]]
  trans <- trans[[1]]
  gc()
  
  
  ####################################################
  ## Prep fertility/birth and birth incidence inputs
  ####################################################
  birth <- prep_birth_cohorts()
  birthprev <- birth[["Sickle cell disorders"]]
  fert <- birth[[1]]
  birth <- NULL
  gc()

  ## loop over disease and location
  for (dis in c("Sickle cell disorders")) {
    for (loc in c(clist$location_name)) {

        cat(paste0(scen,dis,loc,"\n")); flush.console()
      
        ## identify which dimensions of transition array need to be passed based on index
        subset_male <- which(trans_index$sex == "male" & trans_index$cause_name==dis & trans_index$location_name==loc & trans_index$year %in% c(base_year:(base_year+proj_yrs)))
        subset_female <- which(trans_index$sex == "female" & trans_index$cause_name==dis & trans_index$location_name==loc & trans_index$year %in% c(base_year:(base_year+proj_yrs)))

        ## run model
        out <- sim_markov_chain_array_sex(x0=copy(init[cause_name==dis & location_name==loc,c("sex","age","Susceptible","PrevUnenrolled","Enrolled","LTFU","DeadDis","DeadNoDis"),with=F]),
                                   pf=trans[,,subset_female],pm=trans[,,subset_male],
                                   n_cycles=proj_yrs,xadd=fert[location_name==loc],state_names=state_names,base_year=base_year,bprev=birthprev[location_name==loc],
                                   b_enroll_b,b_enroll_d[scen])
          
        ## compile results for males
        males <- rbindlist(lapply(1:dim(out[[1]])[3],FUN=function(x){
          tmp <- data.table(as.data.frame(out[[1]][,,x]))
          tmp[,year:=base_year+x-1]
          return(tmp)
        }))
  
        ## compile results for females
        females <- rbindlist(lapply(1:dim(out[[2]])[3],FUN=function(x){
          tmp <- data.table(as.data.frame(out[[2]][,,x]))
          tmp[,year:=base_year+x-1]
          return(tmp)
        }))
        
        ## add some characteristics for collapsing later
        res[[paste0(scen,loc)]] <- rbind(males,females)
        res[[paste0(scen,loc)]][,sex:=c(rep("male",nrow(res[[paste0(scen,loc)]])/2),rep("female",nrow(res[[paste0(scen,loc)]])/2))]
        res[[paste0(scen,loc)]][,age:=rep(0:95,nrow(res[[paste0(scen,loc)]])/96)]
        res[[paste0(scen,loc)]][,location_name:=loc]
        res[[paste0(scen,loc)]][,scenario:=scen]

    }
  }
}

res <- rbindlist(res)

## relabel scenarios
res[,scenario:=scenario-1]
res[,scenario:=factor(scenario,levels=c(0:5),labels=c("Baseline","1","2","3","4","5"))]

## summarize/rename results from states
res <- res[,list(S=sum(Susceptible),P=sum(PrevUnenrolled),E=sum(Enrolled),LTFU=sum(LTFU),DeadDis=sum(DeadDis),DeadNoDis=sum(DeadNoDis)),by=c("year","location_name","age","scenario")]

##########################
## create primary results
##########################
tmp <- copy(res[scenario %in% c("Baseline","1","2","3","4","5")])
tmp[,prev_num:=P+E+LTFU] ## prevalence of SCD
tmp[,pop:=S+P+E+LTFU] ## total population alive
tmp <- tmp[,list(prev_num=sum(prev_num),pop=sum(pop),DeadDis=sum(DeadDis)),by=c("year","location_name","scenario")]
tmp[,dead_yr:=c(0,diff(DeadDis))] ## rather than total deaths, how many deaths from SCD were there in the year?
tmp[year==2022,dead_yr:=0]
tmp[,mx:=dead_yr/pop] ## mortality rate
tmp[,prev:=prev_num/pop] ## prevalence %
tmp <- dcast.data.table(tmp,year+location_name~scenario,value.var=c("mx","pop","prev"))

## for deaths averted and person-years gained, there is the complication of which scenario's population to use
## i.e. since there are more births if more ppl survive, should those also count as additional survival? Probably not.
## So, instead, calculate difference in deaths/pop and prevalence %, apply those to the baseline scenario's population
## this method has been used elsewhere
tmp <- melt(tmp,id.vars=c("year","location_name","prev_Baseline","pop_Baseline","mx_Baseline"))
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
tmp[,variable:=as.character(variable)]
tmp[,scenario:=substrRight(variable,1)]
tmp[,variable:=substr(variable,1,nchar(variable)-2)]
tmp <- dcast.data.table(tmp,year+location_name+prev_Baseline+pop_Baseline+mx_Baseline+scenario~variable,value.var="value")
tmp[,deaths_averted:=pop_Baseline*(mx_Baseline-mx)]
tmp[,pys_gained:=pop_Baseline*(prev-prev_Baseline)]

## create table output for paper
tabout <- copy(tmp[scenario %in% c("1","2","3"),list(deaths_averted=sum(deaths_averted),pys_gained=sum(pys_gained),pop_Baseline=sum(pop_Baseline)),by=c("location_name","scenario")])
tabout[,death_rate_averted:=deaths_averted/pop_Baseline*100000]
tabout[,pys_rate_gained:=pys_gained/pop_Baseline*100000]
tabout[,deaths_averted:=deaths_averted/1000]
tabout[,pys_gained:=pys_gained/1000]

tabout <- dcast.data.table(tabout,location_name~scenario,value.var = c("deaths_averted","pys_gained","death_rate_averted","pys_rate_gained"))
write.csv(tabout,"/tmp/pen_plus/tables/results_table.csv",row.names=F)



##########################
## sensitivity analysis
#########################
sens <- copy(tmp[scenario %in% c("1","2","3","4","5"),list(deaths_averted=sum(deaths_averted),pys_gained=sum(pys_gained),pop_Baseline=sum(pop_Baseline)),by=c("scenario")])
sens[,deaths_averted:=deaths_averted/1000]
sens[,pys_gained:=pys_gained/1000]

## make tornado plot (single variable)

# width of columns in plot (value between 0 and 1)
width <- 0.95

## prep data to plot in ggplot using rectangles
toplot <- data.table(data.frame(lower=sens[scenario=="5"]$deaths_averted,upper=sens[scenario=="4"]$deaths_averted))
toplot[,variable:=c("Deaths Averted (thousands)")]
toplot <- melt(toplot,id.vars=c("variable"),variable.name="bound")
toplot[,ymin:=pmin(value,sens[scenario=="1"]$deaths_averted)]
toplot[,ymax:=pmax(value,sens[scenario=="1"]$deaths_averted)]
toplot[,xmin:=1-width/2]
toplot[,xmax:=1+width/2]

gg <- ggplot() + 
  geom_rect(data = toplot, 
            aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=bound)) +
  theme_bw() + 
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = sens[scenario=="1"]$deaths_averted) +
  scale_x_continuous(breaks = c(1),
                     labels = c("")) +
  scale_y_continuous("Deaths Averted (thousands)",limits=c(200,1200)) +
  scale_fill_discrete(labels=c("50% Effect Size","90% Effect Size")) +
  coord_flip()
print(gg)




stop()

## data for diagnostic/calibration plots (This was initially run with certain parameters to check calibration)
## RR = 1 to just have derived GBD parameters to check calibration
dat <- copy(res[location_name=="Nigeria"])
dat[,prev_num:=P+E+LTFU]
dat[,pop:=S+P+E+LTFU]
dat[,prev_pct:=prev_num/pop]
dat[,year:=factor(year)]

tmp <- copy(dat[,c("year","age","S","P","E"),with=F])
tmp <- melt(tmp,id.vars=c("year","age"))
tmp <- tmp[,list(value=sum(value)),by=c("year","variable")]
tmp[,year:=as.numeric(year)]

gg <- ggplot(tmp,aes(x=year,y=value)) + geom_line() + facet_wrap(~variable,scales="free")

print(gg)

## prevalence by age over time
pdf("/tmp/pen_plus/figures/NGA_matrix_prevalence.pdf",width=7,height=4)
gg <- ggplot(data=dat[year %in% c(2022,2027,2032,2037,2042)],aes(x=age,y=prev_pct*100,group=year,color=year)) + geom_line(size=1.1) + theme_bw() + 
  ylab("Prevalence of SCD (%)") + xlab("Age (years)") + scale_color_discrete("Year")
print(gg)
dev.off()


## diagnostic, see how far off pop is from projections
pdf("/tmp/pen_plus/figures/NGA_matrix_pop.pdf",width=6,height=4)
gg <- ggplot(data=dat[year %in% c(2022,2042)],aes(x=age,y=pop/1000000,group=year,color=year)) + geom_line(size=1.1) + theme_bw() + 
  ylab("Population (millions)") + xlab("Age (yrs)") + scale_color_discrete("Year")
print(gg)
dev.off()
sum(dat[year==2042]$pop)/1000000 ## 352 million with initial params, UN WPP says 343 million
## differences: GBD pop baseline, migration, cause-deleted mortality vs. all-cause in susceptible pop


# tot <- copy(dat)
# tot <- tot[,list(S=sum(S),P=sum(P),E=sum(E),LTFU=sum(LTFU),DeadDis=sum(DeadDis),DeadNoDis=sum(DeadNoDis),prev_num=sum(prev_num),pop=sum(pop)),by=c("year","location_name")]
# tot <- melt(tot,id.vars=c("year","location_name","pop"))

# pdf("/tmp/pen_plus/figures/NGA_matrix_cases_care.pdf",width=6,height=4)
# 
# gg <- ggplot(data=tot[variable %in% c("P","E")],aes(x=year,y=value/1000000,fill=variable)) + geom_bar(stat="identity",position="stack") + theme_bw() + 
#   ylab("People (millions)") + xlab("Year") + scale_fill_discrete("",labels=c("Unenrolled","Enrolled")) + 
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# print(gg)
# 
# dev.off()


