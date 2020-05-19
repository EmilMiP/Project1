library(dplyr)
library(truncnorm)
library(mvtnorm)

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



omniscience = function(input, h2, thr, nsib = 0, NA_frac, init_samp = 1e7, account_for_sex = FALSE, downsample_frac = 1) { #= 1/(1 + 2*nrow(thr))

  ###assuming names are fixed for now, but might need to make them variable at a later point  
  
#### Simulate multivariate liabilities ####
if (nsib <= 0) {
  cov <- diag(c(h2, 1 - h2, 1, 1))
  colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", "father_full", "mother_full")
  cov[2 + 1:2, 1] <- cov[1, 2 + 1:2] <- h2 / 2 
  cov
#  round(100 * cov(cbind(child_gen, child_env = child_full - child_gen, 
#                        father_full, mother_full)), 2)
} else {### This handling of the siblings might be a major assumption
###  cov <- diag(c(h2, 1 - h2, rep(1 - h2, nsib), 1, 1))
###  colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", if ( nsib > 0) paste("sib",1:nsib, "_full", sep = ""),"father_full", "mother_full")
###  #fixing cov for gen child:
###  cov[2 + nsib + 1:2, 1] <- cov[1, 2 + nsib + 1:2] <- h2 / sqrt(2)  # verify this but seems okay
###  #fixing gen cov for siblings full and gen child:
###  cov[2 + 1:nsib, 1] <- cov[1, 2 + 1:nsib] <- h2  # verify this
###  # sorting out the correlations of the siblings
###  cov[2 + 1:nsib, 2 + 1:nsib] <- cov[2 + 1:nsib, 2 + 1:nsib] + h2 # verify this 
###  #sorting out the correlation between siblings and parents:
###  cov[ncol(cov) - 0:1, ncol(cov) - 1 - 1:nsib] <- cov[ncol(cov) - 1 - 1:nsib, ncol(cov) - 0:1] <- h2 / sqrt(2) # verify this 
###  print(cov)
####  round(100 * cov(cbind(child_gen, 
####                        child_env = child_full - child_gen, 
####                        sib_full = sib_full,
####                        father_full, mother_full)), 2)
}


  if (is.null(dim(thr))) {
    thr_len = length(thr)
  } else {
    thr_len = nrow(thr)
  }

nb_var <- ncol(cov) - 1
if (!account_for_sex) {
  all_config <- do.call(expand.grid, c(rep(list(c(NA,1:(thr_len*2))), nb_var)))
  names(all_config) <- c("child_status", if (nsib > 0) paste("sib",1:nsib, "_status", sep = ""), "father_status", "mother_status")
  all_config$string <- as.factor(do.call(paste, all_config + 0L))
} else {
  all_config <- do.call(expand.grid, c(rep(list(c(NA,1:(thr_len*2))), nb_var), rep(list(1:2), 1 + nsib)))
  names(all_config) <- c("child_status", if (nsib > 0) paste("sib",1:nsib, "_status", sep = ""), "father_status", "mother_status", "child_sex", if (nsib > 0) paste("sib",1:nsib, "_sex", sep = ""))
  all_config$string <- as.factor(do.call(paste, all_config + 0L))
}
all_config = as_tibble(all_config)
cat("Maximum number of configurations:", nrow(all_config), "\n")

#need to make sure that the ordering is the same. otherwise we might compare the wrong columns.
config_cols = sort(str_subset(colnames(input), "status|sex"))
observed_configs = unique(do.call(paste, input[, config_cols]))
cat("Observed number of configurations:", length(observed_configs), "\n")
all_conf_cols = sort(colnames(all_config)[-ncol(all_config)])
all_conf_string = do.call(paste, all_config[,all_conf_cols])
all_config = all_config[all_conf_string %in% observed_configs,]

cat("simulating the initial data. \n")
set.seed(1)
simu_liab <- mvtnorm::rmvnorm(init_samp, sigma = cov)


for (i in 2:ncol(cov)) { #starting at 2 to not filter on os twice, i.e. environment and genetic filtering
  simu_liab[base::sample(1:nrow(simu_liab), nrow(simu_liab) * NA_frac),i] <- NA
}




#seq_n <- seq_len(nrow(simu_liab))

#are we accounting for sex or not ?
if (!account_for_sex) {
  sex = rep(1, ncol(simu_liab) - 1)
} else {
  sex = list(
    child_sex  = sample(1:2, size = nrow(simu_liab), replace = TRUE),
    father_sex = 1,
    mother_sex = 2
  )
  if (nsib > 0) {
    for (ii in 1:nsib) {
      sex[[paste("sib", ii, "_sex", sep = "")]] = sample(1:2, size = nrow(simu_liab), replace = TRUE)
    }
  }
}

#get full liab for child:
simu_liab[,2] = simu_liab[,1] + simu_liab[,2]

all_split <- lapply(2:(3 + 1 + nsib), function(k) {
  print(k)
  #flag for account for sex is fixed earlier. 
  sex_val = sex[[k - 1]]
  if (length(sex_val) > 1) { #sex is different or the same, but it is there.
    case_grps = numeric(nrow(simu_liab)) - 9
    ctrl_grps = numeric(nrow(simu_liab)) - 9
    for (sex_k in 1:2) {
      sex_index = sex_val == sex_k
      case_grps[sex_index] = rowSums(outer(simu_liab[sex_index,k], thr[, sex_k], `>`)) + thr_len
      ctrl_grps[sex_index] = rowSums(outer(simu_liab[sex_index,k], thr[, sex_k], `>`)) + 1
    }
  } else {
    case_grps = rowSums(outer(simu_liab[,k], thr[, sex_val], `>`)) + thr_len
    ctrl_grps = rowSums(outer(simu_liab[,k], thr[, sex_val], `>`)) + 1
  }
  indiv_split = list()
  indiv_split[['0']] = which(is.na(case_grps))
  for (kk in 1:(thr_len)) {
    indiv_split[[paste(kk)]] = which(ctrl_grps <= kk)
  }
  for (kk in (thr_len + 1):(2 * thr_len)) {
    indiv_split[[paste(kk)]] = which(case_grps >= kk)
  }
  indiv_split
})
for (i in 1:(3 + nsib)) {
  all_split[[i]][["sex"]] = sex[[i]]
}


cat("Assigning simulated family to a group: \n")

no_conf = nrow(all_config)
cat("Number of configurations:", no_conf, "\n")

update_point = round(no_conf/20)

#updating all_config to have 0 in stead of NA values
all_config[is.na(all_config)] <- 0
#this step becomes slow if the number of groups is high (!)
all_config$obs <- lapply(1:nrow(all_config), function(i) {
  if (i %% update_point == 0) cat("~", round(i/no_conf * 100), "%\n")
  tmp <- all_config[i,grep("status|sex", colnames(all_config))]
  all_ind <- lapply(1:(3 + nsib), function(k) {
    all_split[[k]][[tmp[[k]] + 1]]
  })
  if (account_for_sex) {
    all_ind[["sex"]] = (1:nrow(simu_liab))[sex[[1]] == tmp$child_sex]
    if (nsib > 0) {
      for (ik in 1:nsib) { ### Still needs to be tested
        sib_sex_name = paste("sib", ik,"_sex", sep = "")
        all_ind[[sib_sex_name]] = (1:nrow(simu_liab))[sex[[ik + 3]] == tmp[[sib_sex_name]]]
      }
    }
  }
  simu_liab[Reduce(intersect, all_ind),1]
})
all_config[["post_mean_liab"]] = sapply(all_config$obs, FUN = mean)
all_config[["se"]] = sapply(all_config$obs, FUN = function(x) sd(x)/sqrt(length(x)))
#all_config[all_config$se > 0.01,]

#### Assign group posterior mean genetic liabilities to individuals ####

child_group <- input %>%
  sample_frac(downsample_frac, weight = ifelse(input$child_status > thr_len, 1e6, 1)) %>%
  left_join(all_config) # %>% print()
return(child_group[,-grep("obs", colnames(child_group))])
}

#omni_dat = omniscience(input, h2, thr, nsib)
#sum(omni_dat$n < 2)







