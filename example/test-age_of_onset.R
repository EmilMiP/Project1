library(dplyr)
library(MASS)
#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 1e4
prev <- 0.05
downsample_frac <- 1 #1 for no downsampling.

set.seed(1) #set seed for more reliable comparison

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

#### Simulate multivariate liabilities ####
covSim <- diag(c(h2, 1 - h2, 1, 1))
colnames(covSim) <- rownames(covSim) <- c("child_gen", "child_env", "father_full", "mother_full")
covSim[3:4, 1] <- covSim[1, 3:4] <- h2 / 2  

#Only need to simulate liabilites once.  We can reuse for different values of K
simu_liab <- mvtnorm::rmvnorm(5e7, sigma = covSim)  

#### Assign case-control status ####
K <- prev * 10:1/10
#K <- prev * 5:1/5
#K <- prev * c(5,1)/5
#K <- c(0.05, 0.01)
#K <- prev
(thr <- qnorm(1 - K))

fStat <- rowSums(outer(fLiab, thr, `>`))
mStat <- rowSums(outer(mLiab, thr, `>`))
cStat  <- rowSums(outer(cLiab, thr, `>`))

nb_var <- ncol(covSim) - 1
all_config <- do.call(expand.grid, rep(list(0:length(thr)), nb_var))
names(all_config) <- c("child_status", "father_status", "mother_status")
all_config$string <- as.factor(do.call(paste, all_config + 0L))

df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_status  = (rowSums(outer(simu_liab[, 1] + simu_liab[, 2], thr, `>`))),
  father_status = (rowSums(outer(simu_liab[, 3], thr, `>`))),
  mother_status = (rowSums(outer(simu_liab[, 4], thr, `>`)))
) %>%
  left_join(all_config) 

group_means <- group_by(df_simu_liab, string, .drop = FALSE) %>%
  summarise(post_mean_liab = mean(child_gen), n = n(),
            se = sd(child_gen) / sqrt(n)) %>%
  ungroup() %>%
  arrange(post_mean_liab) 

#### Assign group posterior mean genetic liabilities to individuals ####
child_group <- tibble(child_gen = cGenLiab, child_status = cStat, father_status = fStat, mother_status = mStat) %>%
  sample_frac(downsample_frac, weight = ifelse(child_status, 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) 

trues = child_group %>% group_by(string) %>% summarise(true_mean_liab = mean(child_gen)) %>% left_join(group_means)

library(ggplot2)
ggplot(child_group) +
#  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen), alpha = 0.2) + 
  geom_point(data = trues, aes(post_mean_liab, true_mean_liab), col = "red") +
  geom_abline(col = "red")

#Squared correlation is a better measurement
c(cor(child_group$post_mean_liab, child_group$child_gen)**2,
  cor((child_group$child_status > 0) + 0L, child_group$child_gen)**2)

### 5% cases + N = 5k
# k = 0.05:           0.4266394 0.3332564
# k = c(0.05, 0.01):  0.4296884 0.3332564
# k = 5:1/100         0.4313874 0.3332564
# k = 10:1/200        0.4319736 0.3332564

### 5% cases + N = 10k
# k = 0.05:           0.4337232 0.3183950
# k = c(0.05, 0.01):  0.4372176 0.3183950
# k = 5:1/100         0.4391433 0.3183950
# k = 10:1/200        0.4390453 0.3183950

### 50% cases + N = 5k
# k = 0.05:           0.7152612 0.6662278
# k = c(0.05, 0.01):  0.7255952 0.6662278
# k = 5:1/100         0.7286037 0.6662278
# k = 10:1/200        0.7305081 0.6662278

### 50% cases + N = 10k
# k = 0.05:           0.7148283 0.6697706
# k = c(0.05, 0.01):  0.7231299 0.6697706
# k = 5:1/100         0.7262446 0.6697706
# k = 10:1/200        0.7278923 0.6697706

