set.seed(42)

# Loading necessary libraries
library(tidyverse)

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

# Some basic constants
h2 <- 0.5
max_prev <- 0.5
sex_eff <- 0.5
all_ages <- seq_len(100)
all_years <- 1920:2019
all_sex <- 1:2

# 1. Simulate cumulative incidence rates

# Can be updated for a non-normal distribution, say using observed simulated distribution.
get_thres <- function(prev) {
  qnorm(prev, lower.tail = FALSE)
}

get_prev <- function(ages, years, sex) {
  prev <- max_prev * 
    (ages - min(all_ages)) / length(all_ages) *
    (years - min(all_years)) / length(all_years)
  ifelse(sex == 2, prev * sex_eff, prev)
}

get_thres_table <- function() {
  expand_grid(age = all_ages, year = all_years, sex = all_sex) %>%
    filter(age < (2020 - year)) %>% # Age is limited by year of birth
    mutate(prev = get_prev(age, year, sex),
           thres = get_thres(prev))
}

thres_table <- get_thres_table()

filter(thres_table, age == 40, year == 1980)

get_cov <- function(h2, n_sib = 0) {
  cov <- matrix(h2 / 2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3, 4] <- cov[4, 3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  cov
}

# 2. Simulate sample
#### Simulating some data (liabilities) ####
simTrioLiab <- function(N, h2) {
  
  cov <- get_cov(h2)
  liab <- MASS::mvrnorm(N, rep(0, ncol(cov)), cov)
  
  tibble(
    child_gen   = liab[, 1],
    child_full  = liab[, 2],
    child_sex   = sample(1:2, N, replace = TRUE),
    child_age   = sample(1:35, N, replace = TRUE),
    child_year  = sample(1981:2010, N, replace = TRUE),
    mother_full = liab[, 3],
    mother_sex  = 1,
    mother_age  = sample(20:80, N, replace = TRUE),
    mother_year = sample(1950:1990, N, replace = TRUE),
    father_full = liab[, 4],
    father_sex  = 2,
    father_age  = sample(20:80, N, replace = TRUE),
    father_year = sample(1950:1990, N, replace = TRUE)
  )
}

trio_liabs <- simTrioLiab(1e6, h2)
# Age is limited by year of birth, etc
trio_liabs <-
  filter(
    trio_liabs,
    child_age  < (2020 - child_year),
    mother_age < (2020 - mother_year),
    father_age < (2020 - father_year),
    child_age < (mother_age - 15),
    child_age < (father_age - 15)
  )
trio_liabs
#Filter impossible combinations

# 3. For each sampled individual identify case status and age-of-onset

get_age_of_onset <- function(l, y, s) {
  ind <- with(thres_table, which(year == y & sex == s & thres <= l))
  if (length(ind) == 0) NA else min(thres_table$age[ind])
}
get_ao <- Vectorize(get_age_of_onset)

get_ea <- function(age, age_of_onset) {
  ifelse(is.na(age_of_onset), age, age_of_onset)
}


sample_table <- trio_liabs %>%
  # add threshold & prevalence for child
  left_join(thres_table,
            by = c(
              "child_age" = "age",
              "child_year" = "year",
              "child_sex" = "sex"
            )) %>%
  mutate(child_status = as.integer(child_full > thres),
         child_age_of_onset = get_ao(child_full, child_year, child_sex), # Not very efficient
         child_eff_age = get_ea(child_age, child_age_of_onset)) %>%
  rename(child_thres = thres, child_prev = prev) %>%
  # add threshold & prevalence for mother
  left_join(thres_table,
            by = c(
              "mother_age" = "age",
              "mother_year" = "year",
              "mother_sex" = "sex"
            )) %>%
  mutate(mother_status = as.integer(mother_full > thres),
         mother_age_of_onset = get_ao(mother_full, mother_year, mother_sex), # Not very efficient
         mother_eff_age = get_ea(mother_age, mother_age_of_onset)) %>%
  rename(mother_thres = thres, mother_prev = prev) %>%
  # add threshold & prevalence for father
  left_join(thres_table,
            by = c(
              "father_age" = "age",
              "father_year" = "year",
              "father_sex" = "sex"
            )) %>%
  mutate(father_status = as.integer(father_full > thres),
         father_age_of_onset = get_ao(father_full, father_year, father_sex), # Not very efficient
         father_eff_age = get_ea(father_age, father_age_of_onset)) %>%
  rename(father_thres = thres, father_prev = prev)

#### Florian's code using the truncated multivariate normal distribution ####
sample_table_cases <-
  filter(sample_table, child_status==1)%>% 
  sample_n(., 1000)
sample_table_controls <-
  filter(sample_table, child_status==0) %>% 
  sample_n(., 2000)
sample_table <- rbind(sample_table_cases,sample_table_controls)


sample_table$post_gen_liab <- NA
sample_table$post_gen_liab2 <- NA
inv_cov <- solve(cov <- get_cov(h2))

max_u_bound = Inf
min_l_bound = -Inf

library(matrixStats)
for (i in 1:nrow(sample_table)) {
  print(i)
  lims <- cbind(upper = rep(max_u_bound, 4), lower = rep(min_l_bound, 4))
  x <- sample_table[i, ]
  if (x$child_status==1){
    lims[2, 1] <- x$child_thres+0.01
    lims[2, 2] <- x$child_thres-0.01
  }
  else{
    lims[2, 1] <- x$child_thres
  }
  if (x$mother_status==1){
    lims[3, 1] <- x$mother_thres+0.01
    lims[3, 2] <- x$mother_thres-0.01
  }
  else{
    lims[3, 1] <- x$child_thres
  }
  
  if (x$father_status==1){
    lims[4, 1] <- x$father_thres+0.01
    lims[4, 2] <- x$father_thres-0.01
  }
  else{
    lims[4, 1] <- x$father_thres
  }
  # lims[3, x$mother_status + 1] <- x$mother_thres
  # lims[4, x$father_status + 1] <- x$father_thres
  liabs <- tmvtnorm::rtmvnorm(2e4, mean = rep(0, 4), H = inv_cov,
                              lower = lims[, "lower"], upper = lims[, "upper"],
                              algorithm = "gibbs", burn.in.samples = 1000)
  
  # liabs <- TruncatedNormal::rtmvnorm(1e3, mu = rep(0, 4), sigma = cov,
  #                                    lb = lims[, "lower"], ub = lims[, "upper"])
  sample_table$post_gen_liab[i] <- colMeans(liabs,na.rm = TRUE)[1]
  sample_table$post_gen_liab2[i] <- colMedians(liabs,na.rm = TRUE)[1]
}

# for (i in 1:nrow(sample_table)) {
#   print(i)
#   lims <- cbind(upper = rep(Inf, 4), lower = rep(-Inf, 4))
#   x <- sample_table[i, ]
#   lims[2, x$child_status  + 1] <- x$child_thres
#   lims[3, x$mother_status + 1] <- x$mother_thres
#   lims[4, x$father_status + 1] <- x$father_thres
#   fixed <- c(FALSE, x$child_status, x$mother_status, x$father_status)
#   # liabs <- tmvtnorm::rtmvnorm(2e4, mean = rep(0, 4), H = inv_cov,
#   #                             lower = lims[, "lower"], upper = lims[, "upper"],
#   #                             algorithm = "gibbs", burn.in.samples = 2000)
#   # liabs <- TruncatedNormal::rtmvnorm(10e3, mu = rep(0, 4), sigma = cov, 
#   #                                    lb = lims[, "lower"], ub = lims[, "upper"])
#   # gen_liabs <- liabs[, 1]
#   gen_liabs <- rtmvnorm.gibbs(50e3, burn_in = 1000, sigma = cov, 
#                               lower = lims[, "lower"], 
#                               upper = lims[, "upper"],
#                               fixed = fixed)
#   # gen_liabs <- gen_liabs[seq(1, length(gen_liabs), by = 20)]  # thinning
#   sample_table$post_gen_liab[i] <- mean(gen_liabs)
#   sample_table$post_gen_liab_se[i] <- sd(gen_liabs) / sqrt(length(gen_liabs))
# }
# summary(sample_table$post_gen_liab_se)

# microbenchmark::microbenchmark(
#   tmvtnorm::rtmvnorm(10e3, mean = rep(0, 4), H = inv_cov, 
#                      lower = lims[, "lower"], upper = lims[, "upper"],
#                      algorithm = "gibbs", burn.in.samples = 1000)[, 1],
#   TruncatedNormal::rtmvnorm(10000, mu = rep(0, 4), sigma = cov, 
#                             lb = lims[, "lower"], ub = lims[, "upper"])[, 1],
#   rtmvnorm.gibbs(10e3, burn_in = 1000, sigma = cov, 
#                  lower = lims[, "lower"], upper = lims[, "upper"]),
#   tmvtnorm::rtmvnorm(10e3, mean = rep(0, 4), H = inv_cov, 
#                      lower = lims[, "lower"], upper = lims[, "upper"])[, 1],
#   times = 10
# )
# 
# Ngibbs <- 1e6
# lapply(list(
#   tmvtnorm::rtmvnorm(Ngibbs, mean = rep(0, 4), H = inv_cov, 
#                      lower = lims[, "lower"], upper = lims[, "upper"],
#                      algorithm = "gibbs", burn.in.samples = 1000)[, 1],
#   TruncatedNormal::rtmvnorm(Ngibbs, mu = rep(0, 4), sigma = cov, 
#                             lb = lims[, "lower"], ub = lims[, "upper"])[, 1],
#   rtmvnorm.gibbs(Ngibbs, burn_in = 1000, sigma = cov, 
#                  lower = lims[, "lower"], upper = lims[, "upper"]),
#   tmvtnorm::rtmvnorm(Ngibbs, mean = rep(0, 4), H = inv_cov, 
#                      lower = lims[, "lower"], upper = lims[, "upper"])[, 1]
# ), summary)

sample_table %>%
  filter(child_status == 1, mother_status == 1, father_status == 0) %>%
  select(child_gen, child_year, child_age, child_prev, child_thres, post_gen_liab) %>%
  print() %>%
  ggplot(aes(child_prev, post_gen_liab)) + 
  geom_point()

ggplot(sample_table) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_gen_liab, child_gen, color = factor(child_status)), alpha = 0.4) + 
  geom_smooth(aes(x=post_gen_liab, y=child_gen)) +
  geom_abline(col = "black", linetype = 2) + 
  theme(legend.position = "top")
with(sample_table, cor(post_gen_liab, child_gen)**2)  
# 0.3345276 with 1e3 -> 0.3405035 with 100e3 / 0.3409902 with 10e3 and 1e3 burn-in

ggplot(sample_table) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_gen_liab, child_gen, color = child_age), alpha = 0.4) + 
  geom_smooth(aes(x=post_gen_liab, y=child_gen)) +
  geom_abline(col = "black", linetype = 2) + 
  theme(legend.position = "top") +
  scale_color_viridis_c()

# standard LT-FH
simu_liab <- mvtnorm::rmvnorm(1e7, sigma = cov)
thr <- median(c(sample_table$father_thres, sample_table$mother_thres))

all_config <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), 3))
names(all_config) <- c("child_status", "father_status", "mother_status")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config

library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_status  = (simu_liab[, 2] > thr),
  father_status = (simu_liab[, 3] > thr),
  mother_status = (simu_liab[, 4] > thr)
) %>%
  left_join(all_config) %>%
  print()

group_means <- group_by(df_simu_liab, string, .drop = FALSE) %>%
  summarise(post_mean_liab = mean(child_gen), n = n(),
            se = sd(child_gen) / sqrt(n)) %>%
  ungroup() %>%
  arrange(post_mean_liab) %>%
  print()

#### Assign group posterior mean genetic liabilities to individuals ####
child_group <- sample_table %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()

ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen, color = factor(child_sex)), alpha = 0.4) + 
  geom_abline(col = "black", linetype = 2) + 
  geom_smooth(aes(x=post_gen_liab, y=child_gen)) +
  theme(legend.position = "top")
with(child_group, cor(post_mean_liab, child_gen)**2)  # 0.3394834


#### Bjarni's original code below ####

############ -- This concludes the input data simulation -- #############
#
# 4. For each genotyped individual identify case-control status identify integral region
#    a) Group individuals by
#         - age (for controls)
#         - age-of-onset (for cases)
#         - year
#         - sex

year_group_thres <- c(1990, 2000)
year_groups <- seq(length(year_group_thres) + 1)
age_group_thres <- c(30, 60)
age_groups <- seq(length(age_group_thres) + 1)
sex_groups <- c(1, 2)

group_table = expand_grid(year_group = year_groups,
                          age_group = age_groups,
                          sex_group = sex_groups)
group_table

map_2_age_groups <- function(v) {
  return(sum(v > age_group_thres) + 1)
}
m2ag <- Vectorize(map_2_age_groups)

map_2_year_groups <- function(v) {
  return(sum(v > year_group_thres) + 1)
}
m2yg <- Vectorize(map_2_year_groups)

sample_table_ext <- sample_table %>%
  mutate(child_age_group = m2ag(child_eff_age)) %>%
  mutate(mother_age_group = m2ag(mother_eff_age)) %>%
  mutate(father_age_group = m2ag(father_eff_age)) %>%
  mutate(child_year_group = m2yg(child_year)) %>%
  mutate(mother_year_group = m2yg(mother_year)) %>%
  mutate(father_year_group = m2yg(father_year))
sample_table_ext

#    b) Determine prevalences in each group
##Get observed prevalences and thresholds for each group
thres_table <- thres_table %>%
  mutate(age_group = m2ag(age)) %>%
  mutate(year_group = m2yg(year)) %>%
  group_by(age_group, year_group, sex) %>%
  mutate(group_prev = mean(prev), ) %>%
  mutate(group_thres = qnorm(1 - group_prev), )

thres_groups <- thres_table %>%
  select(sex, age_group, year_group, group_prev, group_thres) %>%
  distinct()

#    c) For each individual identify integral thresholds based on group
#         - Note that a case in age-of-onset group 2 will be integrated between 2 and 1.
# Loop over all individuals

get_int_interval <- function(status, ag, yg, sg) {
  group_thres <-
    filter(thres_groups, age_group == ag, year_group == yg, sex == sg)[[1, "group_thres"]]
  if (status == 0) {
    lbound <- -Inf
    ubound <- group_thres
  }
  else{
    lbound <- group_thres
    cond_age_groups = filter(thres_groups, year_group == yg, sex == sg)["age_group"]
    if (ag == max(cond_age_groups)) {
      ubound <- Inf
    }
    else{
      ubound <-
        filter(thres_groups, age_group == ag + 1, year_group == yg, sex == sg)[[1, "group_thres"]]
    }
  }
  return (list("lbound" = lbound, "ubound" = ubound))
}


clbound = c()
cubound = c()
mlbound = c()
mubound = c()
flbound = c()
fubound = c()
for (i in 1:nrow(sample_table_ext)) {
  #child_thresholds
  cstatus = sample_table_ext[[i, "child_status"]]
  cag = sample_table_ext[[i, "child_age_group"]]
  cyg = sample_table_ext[[i, "child_year_group"]]
  csg = sample_table_ext[[i, "child_sex"]]
  l = get_int_interval(cstatus, cag, cyg, csg)
  clbound[i] = l$lbound
  cubound[i] = l$ubound
  
  mstatus = sample_table_ext[[i, "mother_status"]]
  mag = sample_table_ext[[i, "mother_age_group"]]
  myg = sample_table_ext[[i, "mother_year_group"]]
  msg = sample_table_ext[[i, "mother_sex"]]
  l = get_int_interval(mstatus, mag, myg, msg)
  mlbound[i] = l$lbound
  mubound[i] = l$ubound
  
  fstatus = sample_table_ext[[i, "father_status"]]
  fag = sample_table_ext[[i, "father_age_group"]]
  fyg = sample_table_ext[[i, "father_year_group"]]
  fsg = sample_table_ext[[i, "father_sex"]]
  l = get_int_interval(fstatus, fag, fyg, fsg)
  flbound[i] = l$lbound
  fubound[i] = l$ubound
}

int_interval_table = tibble(clbound, cubound, mlbound, mubound, flbound, fubound)

sample_table_ext <- cbind(sample_table_ext, int_interval_table)


# 5. Simulate liabilities and estimate posterior means
#    a) Simulate liabilities
#    b) Identify possible combinations of thresholds (presumably thousands) and create an ordered list
#         - e.g. ordered by first threshold in each dimension
#    c) For each simulated liability identify groups that it is in and add genetic liability
#    d) Calculate average liabilities in each group
#    e) If not enough simulations, then repeat 4 (or use importance sampling).

#### Simulate multivariate liabilities ####
covSim <- diag(c(h2, 1 - h2, 1, 1))
colnames(covSim) <-
  rownames(covSim) <-
  c("child_gen", "child_env", "father_full", "mother_full")
covSim[3:4, 1] <- covSim[1, 3:4] <- h2 / 2

#Only need to simulate liabilites once.  We can reuse for different values of K
N_sim = 5e3
simu_liab <- mvtnorm::rmvnorm(N_sim, sigma = covSim)
config_table = distinct(int_interval_table)
config_table <- config_table %>%
  add_column(liab_sum = rep(0, nrow(config_table))) %>%
  add_column(num_hits = rep(0, nrow(config_table)))
config_table

for (i in 1:nrow(simu_liab)) {
  cl <- simu_liab[i, 1] + simu_liab[i, 2]
  ml <- simu_liab[i, 3]
  fl <- simu_liab[i, 4]
  for (j in 1:nrow(config_table)) {
    in_reg <-
      cl > config_table[[j, 1]] &
      cl < config_table[[j, 2]] &
      ml >= config_table[[j, 3]] &
      ml < config_table[[j, 4]] &
      fl >= config_table[[j, 5]] & fl < config_table[[j, 6]]
    if (in_reg) {
      config_table[[j, 7]] = config_table[[j, 7]] + simu_liab[i, 1]
      config_table[[j, 8]] = config_table[[j, 8]] + 1
    }
  }
  print(config_table[, "num_hits"])
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
