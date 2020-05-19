### implements the variable threshold for both cases and controls, rather than cases OR(!) controls.

#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 5e4
prev <- c(0.08, 0.02) * 4#.1
downsample_frac <- .4 #1 for no downsampling.
nsib = 0
account_for_sex = FALSE
NA_frac = .2
K <- matrix(prev,ncol = 2)
(thr <- qnorm(1 - K))

#par(mfrow = c(2,3))


# Assisting Functions: ----------------------------------------------------
assign_status = function(input, thr) {
  len_thr = length(thr)
  #check relation to thresholds
  input_status = rowSums(outer(input, thr, `>`)) + len_thr
  input_ctrl = which(input_status == len_thr)
  input_status[input_ctrl] = 1
  return(input_status)
}
assign_status_sex = function(input, thr, input_sex) {
  N_input = length(input) #input_sex = sample(1:2, N_input, replace = TRUE)
  input_status = numeric(N_input)
  input_status[input_sex == 1] = assign_status(input = input[input_sex == 1], thr[,1])
  input_status[input_sex == 2] = assign_status(input[input_sex == 2], thr[,2])
  return(input_status)
}





#### Simulate multivariate liabilities ####
if (nsib <= 0) {
  cov <- diag(c(h2, 1 - h2, 1, 1))
  colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", "father_full", "mother_full")
  cov[2 + 1:2, 1] <- cov[1, 2 + 1:2] <- h2 / 2
  cov
} else {### This handling of the siblings might be a major assumption
  cov <- diag(c(h2, 1 - h2, rep(1 - h2 / 2, nsib), 1, 1))
  colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", if ( nsib > 0) paste("sib",1:nsib, "_full", sep = ""),"father_full", "mother_full")
  #fixing cov for gen child:
  cov[2 + nsib + 1:2, 1] <- cov[1, 2 + nsib + 1:2] <- h2 / 2
  #fixing gen cov for siblings full and gen child:
  cov[2 + 1:nsib, 1] <- cov[1, 2 + 1:nsib] <- h2 / 2
  # sorting out the correlations of the siblings
  cov[2 + 1:nsib, 2 + 1:nsib] <- cov[2 + 1:nsib, 2 + 1:nsib] + h2 / 2 
  #sorting out the correlation between siblings and parents:
  cov[ncol(cov) - 0:1, ncol(cov) - 1 - 1:nsib] <- cov[ncol(cov) - 1 - 1:nsib, ncol(cov) - 0:1] <- h2 / 2 
  cov
}


set.seed(1)
family = rmvnorm(n = N, mean = rep(0, 4 + nsib), sigma = cov)
colnames(family) <- colnames(cov)

input <- tibble(
  ID           = 1:N,
  child_gen    = family[,1],
  child_full   = family[,1] + family[,2],
  child_sex    = sample(1:2, N, replace = TRUE),
  father_full  = family[,ncol(family) - 1],
  mother_full  = family[, ncol(family)]
)


if ( nsib > 0) { ### This handling of the siblings might be a major assumption, i.e. sharing genetic liability with the siblings
  for (i in 1:nsib) {
    input[[paste("sib", i ,"_full", sep = "")]] = family[,2 + i]
    input[[paste("sib", i, "_sex", sep = "")]] = sample(1:2, N, replace = TRUE)
  }
}




#### Assign case-control status ####

input[["father_status"]] = assign_status(input$father_full, thr[,1])
input[["mother_status"]] = assign_status(input$mother_full, thr[,2])
input[["child_status"]]  = assign_status_sex(input$child_full, thr, input$child_sex)
if ( nsib > 0) {
  for (i in 1:nsib) {
    input[[paste("sib", i, "_status", sep = "")]] = input[[4 + 2 * i]] > thr[input[[5 + 2 * i]]]
  }
}

#did not work, "cannot allocate vector of size 372.g Gb"
#child_status  <- rowSums(outer(child_full, thr[,child_sex], `>`))
if (!account_for_sex) thr <- matrix(qnorm(1 - apply(K, MARGIN = 1, mean)), ncol = 2, nrow = nrow(K)) #if all thresholds are the same, then some groups will not get any assigned to them
#if (!account_for_age) {
#  thr <- thr[1,]
#}


nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, c(rep(list(c(NA,1:(nrow(K)*2))), nb_var), list(1:2))) ### change
names(all_config) <- c("child_status", "father_status", "mother_status", "child_sex")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config


set.seed(1)
simu_liab <- mvtnorm::rmvnorm(2e6, sigma = cov)


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

child_group <- input %>%
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



#prev = .1

res = matrix(NA, ncol = 3, nrow = 4)

Klist = list(
  matrix(prev, ncol = 2),
  outer(c(5,1)/5, prev, "*"),
  outer(5:1/5, prev, "*")
 # outer(10:1/10, prev, "*")
)


account_for_sex = TRUE

for (i in seq_along(Klist)) {
  (thr <- qnorm(1 - Klist[[i]]))
  
  #  len_thr = length(thr)
  #### Assign case-control status ####
  input[["father_status"]] = assign_status(input$father_full, thr[,1])
  input[["mother_status"]] = assign_status(input$mother_full, thr[,2])
  input[["child_status"]]  = assign_status_sex(input$child_full, thr, input$child_sex)
 
  omni_dat = omniscience(input = input, thr = thr, h2 = 0.5, init_samp = 2e6, account_for_sex = account_for_sex, NA_frac = 1/(1 + length(thr)*2))
  colnames(omni_dat)[11] = "omni"
  ph = child_group %>% left_join(omni_dat, by = c("ID", "child_gen"))
  res[i,] = c(cor(ph$post_mean_liab, ph$child_gen),
              cor(ph$child_status.x, ph$child_gen),
              cor(ph$omni, ph$child_gen))
  p <- ggplot(omni_dat) +
    bigstatsr::theme_bigstatsr() + 
    geom_point(aes(omni, child_gen, 
                   color = as.factor(child_sex)), alpha = 0.3) + 
    #  geom_point(data = trues, aes(post_mean_liab, true_mean_liab)) +
    geom_abline(col = "black") + 
    theme(legend.position = "top")
  print(p)
}
res
res[,3]/res[,1]


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


