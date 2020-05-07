set.seed(42)

#Loading necessary libraries
library(MASS)
library(tidyverse)
library(dplyr)

# 1. Simulate cumulative incidence rates
#    a) By sex, age, and year
#
# 2. Simulate sample
#    a) Simulate liabilities
#    b) Simulate age, sex, year
#
# 3. For each sampled individual identify case status and age-of-onset
#
############ -- This concludes the input data simulation -- #############
# 
# 4. For each genotyped individual identify case-control status identify integral region
#    a) Group individuals by 
#         - age (for controls) 
#         - age-of-onset (for cases)
#         - year
#         - sex
#    b) Determine prevalences in each group
#    c) For each individual identify integral thresholds based on group
#         - Note that a case in age-of-onset group 2 will be integrated between 2 and 1.
#
# 5. Simulate liabilities and estimate posterior means
#    a) Simulate liabilities
#    b) Identify possible combinations of thresholds (presumably thousands) and create an ordered list 
#         - e.g. ordered by first threshold in each dimension
#    c) For each simulated liability identify groups that it is in and add genetic liability
#    d) Calculate average liabilities in each group
#    e) If not enough simulations, then repeat 4 (or use importance sampling).
#
# 6. For each sampled individual, look up the relevant posterior mean.

#Some basic constants
h2 <- 0.5
max_prev <- 0.5
sex_eff <- 0.5
all_ages <- seq(100)
all_years <- c(1920:2019)
all_sex <- c(1,2)

# 1. Simulate cumulative incidence rates

#Can be updated for a non-normal distribution, say using observed simulated distribution.
get_thres <- function(prev){
  thr <- qnorm(1 - prev)
  return(thr)
}
get_prev <- function(ages,years,sex){
  prev = ((ages-min(all_ages))/length(all_ages))*((years-min(all_years))/length(all_years))*max_prev
  prev[sex==2] = prev[sex==2]*sex_eff
  return(prev)
}

get_thres_table <- function(){
  d <-  expand_grid(age=all_ages, year=all_years, sex=all_sex)
  d <- filter(d,age<2020-year) # Age is limited by year of birth
  prev_list = get_prev(d['age'],d['year'],d['sex'])
  d <- cbind(d, prev_list)
  names(d)[4] <- 'prev'
  d <- cbind(d, get_thres(d[,4]))
  names(d)[5] <- 'thres'
  return(d)
}

thres_table <- get_thres_table()

filter(thres_table,age=40,year=1980)

# 2. Simulate sample
#### Simulating some data (liabilities) ####
simTrioLiab = function(N, h2) {
  gCov = diag(h2, 3)
  gCov[3,1] <- gCov[1,3] <- 0.5*h2
  gCov[3,2] <- gCov[2,3] <- 0.5*h2
  gLiab = mvrnorm(N, rep(0, 3), gCov)
  
  eCov = diag(1-h2, 3)
  eLiab = mvrnorm(N, rep(0, 3), eCov)
  
  fGenLiab <- gLiab[,1]
  fEnvLiab <- eLiab[,1]
  mGenLiab <- gLiab[,2]
  mEnvLiab <- eLiab[,2]
  cGenLiab <- gLiab[,3]
  cEnvLiab <- eLiab[,3]
  
  fLiab <- fGenLiab + fEnvLiab 
  mLiab <- mGenLiab + mEnvLiab 
  cLiab <- cGenLiab + cEnvLiab 

  c_sex <- sample(1:2, N, replace = TRUE)
  c_age <- sample(1:35, N, replace = TRUE)
  c_year <- sample(1981:2010, N, replace = TRUE)
  m_sex <- rep(1,N)
  m_age <- sample(20:80, N, replace = TRUE)
  m_year <- sample(1950:1990, N, replace = TRUE)
  f_sex <- rep(2,N)
  f_age <- sample(20:80, N, replace = TRUE)
  f_year <- sample(1950:1990, N, replace = TRUE)
  
  trio_liabs = tibble(iid = seq(N), child_gen = cGenLiab, child_full = cLiab, child_sex = c_sex, child_age = c_age, child_year=c_year,
                      mother_gen = mGenLiab, mother_full=mLiab, mother_sex = m_sex, mother_age = m_age, mother_year=m_year,
                      father_gen = fGenLiab, father_full=fLiab, father_sex = f_sex, father_age = f_age, father_year=f_year)
  return(trio_liabs)
}

trio_liabs <- simTrioLiab(1e5,h2)
# Age is limited by year of birth, etc
trio_liabs <- filter(trio_liabs,child_age<2020-child_year, mother_age<2020-mother_year, father_age<2020-father_year, child_age<mother_age-15, child_age<father_age-15) 
trio_liabs
#Filter impossible combinations

# 3. For each sampled individual identify case status and age-of-onset

get_age_of_onset <- function(l,y,s){
  t = filter(thres_table,year==y,sex==s,thres<=l)
  if (nrow(t)==0){
    return (NA)
  }
  else{
    return(min(t['age']))
  }
}
get_ao <- Vectorize(get_age_of_onset)

get_eff_age <- function(age,age_of_onset){
  if (is.na(age_of_onset)){
    return (age)
  }
  else{
    return(age_of_onset)
  }
}
get_ea <- Vectorize(get_eff_age)

child_info <- trio_liabs %>% 
  left_join(thres_table,by = c("child_age" = "age", 'child_year'='year', 'child_sex'='sex')) %>% 
  select(-starts_with('father'))%>% 
  select(-starts_with('mother'))%>%
  mutate(child_status=as.integer(child_full>thres))%>%
  mutate(child_age_of_onset=get_ao(child_full,child_year,child_sex))%>% #Not very efficient
  mutate(child_eff_age=get_ea(child_age,child_age_of_onset))%>%
  rename(child_thres=thres,child_prev=prev)

t <- filter(child_info,child_status==0)
print(t)

mother_info <- trio_liabs %>% 
  left_join(thres_table,by = c("mother_age" = "age", 'mother_year'='year', 'mother_sex'='sex')) %>% 
  select(-starts_with('father'))%>% 
  select(-starts_with('child'))%>% 
  mutate(mother_status=as.integer(mother_full>thres))%>% 
  mutate(mother_age_of_onset=get_ao(mother_full,mother_year,mother_sex))%>% #Not very efficient
  mutate(mother_eff_age=get_ea(mother_age,mother_age_of_onset))%>%
  rename(mother_thres=thres,mother_prev=prev)

filter(mother_info,mother_status==1)

father_info <- trio_liabs %>% 
  left_join(thres_table,by = c("father_age" = "age", 'father_year'='year', 'father_sex'='sex')) %>% 
  select(-starts_with('mother'))%>% 
  select(-starts_with('child'))%>%
  mutate(father_status=as.integer(father_full>thres))%>%
  mutate(father_age_of_onset=get_ao(father_full,father_year,father_sex))%>% #Not very efficient
  mutate(father_eff_age=get_ea(father_age,father_age_of_onset))%>%
  rename(father_thres=thres,father_prev=prev)
father_info
filter(father_info,father_status==1)

#Now combine into a large table summarizing all information simulated
sample_table = inner_join(child_info,mother_info,by='iid')
sample_table = inner_join(sample_table,father_info,by='iid')

############ -- This concludes the input data simulation -- #############
# 
# 4. For each genotyped individual identify case-control status identify integral region
#    a) Group individuals by 
#         - age (for controls) 
#         - age-of-onset (for cases)
#         - year
#         - sex

year_group_thres <- c(1990,2000)
year_groups <- seq(length(year_group_thres)+1)
age_group_thres <- c(30,60)
age_groups <- seq(length(age_group_thres)+1)
sex_groups <- c(1,2)

group_table = expand_grid(year_group=year_groups,age_group=age_groups,sex_group=sex_groups)
group_table

map_2_age_groups <- function(v){
  return(sum(v>age_group_thres)+1)
}
m2ag <- Vectorize(map_2_age_groups)

map_2_year_groups <- function(v){
  return(sum(v>year_group_thres)+1)
}
m2yg <- Vectorize(map_2_year_groups)

sample_table_ext <-sample_table %>% 
  mutate(child_age_group=m2ag(child_eff_age))%>%
  mutate(mother_age_group=m2ag(mother_eff_age))%>%
  mutate(father_age_group=m2ag(father_eff_age))%>%
  mutate(child_year_group=m2yg(child_year))%>%
  mutate(mother_year_group=m2yg(mother_year))%>%
  mutate(father_year_group=m2yg(father_year))
sample_table_ext

#    b) Determine prevalences in each group
##Get observed prevalences and thresholds for each group
thres_table <- thres_table %>% 
  mutate(age_group=m2ag(age))%>%
  mutate(year_group=m2yg(year))%>%
  group_by(age_group,year_group,sex)%>%
  mutate(
    group_prev = mean(prev),
  ) %>%
  mutate(
    group_thres = qnorm(1 - group_prev),
  )

thres_groups <- thres_table %>% 
  select(sex,age_group,year_group,group_prev,group_thres)%>%
  distinct()

#    c) For each individual identify integral thresholds based on group
#         - Note that a case in age-of-onset group 2 will be integrated between 2 and 1.
# Loop over all individuals

get_int_interval <- function(status,ag,yg,sg){
  group_thres <- filter(thres_groups,age_group==ag,year_group==yg,sex==sg)[[1,'group_thres']]
  if(status==0){
    lbound <- -Inf
    ubound <-group_thres
  }
  else{
    lbound <-group_thres
    cond_age_groups = filter(thres_groups,year_group==yg,sex==sg)['age_group']
    if(ag==max(cond_age_groups)){
      ubound <- Inf
    }
    else{
      ubound <- filter(thres_groups,age_group==ag+1,year_group==yg,sex==sg)[[1,'group_thres']]
    }
  }
  return (list('lbound'=lbound,'ubound'=ubound))
}


clbound = c()
cubound = c()
mlbound = c()
mubound = c()
flbound = c()
fubound = c()
for (i in 1:nrow(sample_table_ext)){
  #child_thresholds
  cstatus = sample_table_ext[[i,'child_status']]
  cag = sample_table_ext[[i,'child_age_group']]
  cyg = sample_table_ext[[i,'child_year_group']]
  csg = sample_table_ext[[i,'child_sex']]
  l = get_int_interval(cstatus,cag,cyg,csg)
  clbound[i] = l$lbound
  cubound[i] = l$ubound

  mstatus = sample_table_ext[[i,'mother_status']]
  mag = sample_table_ext[[i,'mother_age_group']]
  myg = sample_table_ext[[i,'mother_year_group']]
  msg = sample_table_ext[[i,'mother_sex']]
  l = get_int_interval(mstatus,mag,myg,msg)
  mlbound[i] = l$lbound
  mubound[i] = l$ubound

  fstatus = sample_table_ext[[i,'father_status']]
  fag = sample_table_ext[[i,'father_age_group']]
  fyg = sample_table_ext[[i,'father_year_group']]
  fsg = sample_table_ext[[i,'father_sex']]
  l = get_int_interval(fstatus,fag,fyg,fsg)
  flbound[i] = l$lbound
  fubound[i] = l$ubound
}

int_interval_table = tibble(clbound,cubound,mlbound,mubound,flbound,fubound)

sample_table_ext<- cbind(sample_table_ext,int_interval_table)


# 5. Simulate liabilities and estimate posterior means
#    a) Simulate liabilities
#    b) Identify possible combinations of thresholds (presumably thousands) and create an ordered list 
#         - e.g. ordered by first threshold in each dimension
#    c) For each simulated liability identify groups that it is in and add genetic liability
#    d) Calculate average liabilities in each group
#    e) If not enough simulations, then repeat 4 (or use importance sampling).

#### Simulate multivariate liabilities ####
covSim <- diag(c(h2, 1 - h2, 1, 1))
colnames(covSim) <- rownames(covSim) <- c("child_gen", "child_env", "father_full", "mother_full")
covSim[3:4, 1] <- covSim[1, 3:4] <- h2 / 2  

#Only need to simulate liabilites once.  We can reuse for different values of K
N_sim = 5e3
simu_liab <- mvtnorm::rmvnorm(N_sim, sigma = covSim)  
config_table = distinct(int_interval_table)
config_table <-config_table %>%
  add_column(liab_sum = rep(0,nrow(config_table)))%>%
  add_column(num_hits = rep(0,nrow(config_table)))
config_table

for (i in 1:nrow(simu_liab)){
  cl <- simu_liab[i,1]+simu_liab[i,2]
  ml <- simu_liab[i,3]
  fl <- simu_liab[i,4]
  for (j in 1:nrow(config_table)){
    in_reg <- cl>config_table[[j,1]] & cl<config_table[[j,2]] & ml>=config_table[[j,3]] & ml<config_table[[j,4]] & fl>=config_table[[j,5]] & fl<config_table[[j,6]]
    if (in_reg){
      config_table[[j,7]] = config_table[[j,7]]+simu_liab[i,1]
      config_table[[j,8]] = config_table[[j,8]]+1
    }
  }
  print (config_table[,'num_hits'])
}


#### Assign group posterior mean genetic liabilities to individuals ####

# Plot

####/w sex + age, 50% cases
#prev <- c(0.08, 0.02) + N = 10k      0.7161538 0.6231014

####/w age, 50% cases
#prev <- c(0.08, 0.02) + N = 10k      0.7271203 0.6362876


#### LTFH
# prev <- 0.05 , N = 10k,  50% cases  0.7098228 0.6697183

#### LTFH 
# prev <- c(0.05, 0.05) * 2, N = 10k, 50% cases  0.6985646 0.6340386

#all, 
#0.7195118 0.6530479 vs 0.6923003 0.6438172


