set.seed(1)
h2 <- 0.5
N <- 5e4
NA_frac = 1/3

#max_prev is between [0,1]
#ages is defined between [0,1]
#years is defined between [0,1]
#sex is 1 or 2
#sex_eff is defined between [0,1]
get_prev <- function(max_prev,ages,years,sex,sex_eff=1){
  prev = ages*years*max_prev
  prev[sex==2]*sex_eff
  return(prev)
}

#Can be updated for a non-normal distribution, say using observed simulated distribution.
get_thres <- function(prev){
    thr <- qnorm(1 - prev)
    return(thr)
}

get_cc_stat <- function(liab, thr){
  cc_status = liab>thr
  return(cc_status)
}


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
  trio_liabs = tibble(child_gen = cGenLiab, child_full = cLiab, mother_gen = mGenLiab, mother_full=mLiab, father_gen = fGenLiab, father_full=fLiab)
  return(trio_liabs)
}

trio_liabs <- simTrioLiab(N,h2)


#### Simulate multivariate liabilities ####
covSim <- diag(c(h2, 1 - h2, 1, 1))
colnames(covSim) <- rownames(covSim) <- c("child_gen", "child_env", "father_full", "mother_full")
covSim[3:4, 1] <- covSim[1, 3:4] <- h2 / 2  

#Only need to simulate liabilites once.  We can reuse for different values of K
simu_liab <- mvtnorm::rmvnorm(2e7, sigma = covSim)  

#### Now simulate age, year of birth, sex, etc. ####

# Sampling age, year and sex
c_sex <- sample(1:2, N, replace = TRUE)
c_age <- sample(1:2, N, replace = TRUE)/2
c_year <- sample(1:2, N, replace = TRUE)/2
m_sex <- rep(1,N)
m_age <- sample(1:2, N, replace = TRUE)/2
m_year <- sample(1:2, N, replace = TRUE)/2
f_sex <- rep(2,N)
f_age <- sample(1:2, N, replace = TRUE)/2
f_year <- sample(1:2, N, replace = TRUE)/2

#Now get prevelances
max_prev = 0.5
c_prev = get_prev(max_prev,c_age,c_year,c_sex,sex_eff=0.5)
m_prev = get_prev(max_prev,m_age,m_year,m_sex,sex_eff=0.5)
f_prev = get_prev(max_prev,f_age,f_year,f_sex,sex_eff=0.5)

#Now get thresholds
c_thres = get_thres(c_prev)
m_thres = get_thres(m_prev)
f_thres = get_thres(f_prev)

#Now get case-control status
child_status <- get_cc_stat(trio_liabs['child_full'],c_thres)
mother_status <- get_cc_stat(trio_liabs['mother_full'],m_thres)
father_status <- get_cc_stat(trio_liabs['father_full'],f_thres)

#Now figure out all possible configurations based on thresholds/prevalences
nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, c(rep(list(c(NA,1:(nrow(K)*2))), nb_var), list(1:2))) ### change
names(all_config) <- c("child_status", "father_status", "mother_status", "child_sex")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config



for (i in 2:ncol(cov)) { #starting at 2 to not filter on os twice, i.e. environment and genetic filtering
  simu_liab[base::sample(1:nrow(simu_liab), nrow(simu_liab) * NA_frac),i] <- NA
}


library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_sex     = sample(1:2, nrow(simu_liab), replace = TRUE),
  child_status  = assign_status_sex(simu_liab[, 1] + simu_liab[, 2], thr, child_sex),
  father_status = assign_status(simu_liab[, 3], thr[,1]),
  mother_status = assign_status(simu_liab[, 4], thr[,2])
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

child_group <- tibble(child_gen, child_status, father_status, mother_status, child_sex) %>%
  sample_frac(downsample_frac, weight = ifelse(child_status > nrow(thr), 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
count(child_group, child_status, child_sex)
count(child_group, child_status, "GWAS" = (child_status > nrow(thr)) + 0L)

trues = child_group %>% group_by(string) %>% summarise(true_mean_liab = mean(child_gen), n_true = n(), true_se = sd(child_gen) / sqrt(n_true)) %>% left_join(group_means)

library(ggplot2)
ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen, 
                 color = as.factor(child_sex)), alpha = 0.3) + 
  #  geom_point(data = trues, aes(post_mean_liab, true_mean_liab)) +
  geom_abline(col = "black") + 
  theme(legend.position = "top")
with(child_group, c(cor(post_mean_liab, child_gen), cor(child_status > nrow(K), child_gen)))

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


