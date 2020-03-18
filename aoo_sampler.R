library(truncnorm)
library(mvtnorm)
library(dplyr)
library(purrr)
library(stringr)
library(data.table)
library(doSNOW)
library(foreach)
library(parallel)


sem <- function(x){
  if (length(x) > 1) {
    return(sd(x)/sqrt(length(x)))
  } else {
    return(Inf)
  }
}

resampler = function(sample_size = 5e5, thres = thres, s = s, h2 = h2, status.cols, cols = cols) {
  T_val_dad = thres[[1]][[status.cols$dad_status + 1]]
  T_val_mom = thres[[1]][[status.cols$mom_status + 1]]
  T_val_child = thres[[1]][[status.cols$child_status + 1]]
  
  #simulate parents lia. given status:
  dad_cond = rtruncnorm(n = sample_size,a = T_val_dad[1], b = T_val_dad[2] , mean = 0, sd = 1)
  mom_cond = rtruncnorm(n = sample_size,a = T_val_mom[1], b = T_val_mom[2] , mean = 0, sd = 1)
  
  #mean of (e_g, e_s) given (e_p1,e_p2)
  conditional_means = 0.5*h2*(dad_cond + mom_cond)
  
  #Construction of covariance matrix of (e_g, e_s) given (e_p1,e_p2):
  cov_matrix = diag(s + 1)
  cov_matrix[1,1] = h2*(1 - 0.5*h2)
  if (s > 0) {
    S = (diag(1 - 0.5*(h2)^2 - 0.5*h2*(1 - h2),s,s)) + (0.5*h2*(1 - h2))
    buffer = matrix(rep(0.5*h2*(1 - h2),s),nrow = s,ncol = 1)
    buffer_T = t(buffer)
    
    cov_matrix[1,2:(s + 1)] = buffer_T
    cov_matrix[2:(s + 1),1] = buffer
    cov_matrix[2:(s + 1),2:(s + 1)] = S
  }
  #Construction of covariance matrix of (e_g,e_s, e_p1, e_p2) given e_o (full lia. of offspring):
  cov_matrix_parents_eg_sibs = (diag(1 - (0.5*h2)^2 - 0.5*h2*(1 - 0.5*h2),s + 3,s + 3)) + 0.5*h2*(1 - 0.5*h2)
  buffer_eg_sibs = matrix(rep(0.5*h2*(1 - h2),s + 2),nrow = s + 2,ncol = 1)
  cov_matrix_parents_eg_sibs[1,1] = h2*(1 - h2)
  cov_matrix_parents_eg_sibs[1,2:(s + 3)] = t(buffer_eg_sibs)
  cov_matrix_parents_eg_sibs[2:(s + 3),1] = buffer_eg_sibs
  cov_matrix_parents_eg_sibs[(s + 2),(s + 3)] = cov_matrix_parents_eg_sibs[(s + 3),(s + 2)]  = -1*(0.5*h2)^2
  
  # simulation of (e_g, e_s) given (e_p1, e_p2) for all configurations of parental stati:
  if (s > 0) {
    offspring_genetic_sibling_eps  = rmvnorm(n = sample_size, sigma = cov_matrix) + cbind(conditional_means, matrix(rep(conditional_means, s), ncol = s))
    sib.col.names = paste("sib", 1:s, sep = "")
    colnames(offspring_genetic_sibling_eps)[2:(1 + s)] = sib.col.names
  } else {
    offspring_genetic_sibling_eps  = rnorm(n = sample_size, mean = conditional_means, cov_matrix)
  }
  #first set of offspring simulations is done here.
  
  #simulate offsprings lia. given status:
  offspring_epsilon_o = rtruncnorm(n = sample_size, a = T_val_child[1], b = T_val_child[2], mean = 0, sd = 1)
  
  #calculating mean and simulating from (e_g, e_s, e_p1, e_p2) given e_o (full lia).
  #it is done in both configurations:
  offspring_genetic_sibling_parents_eo_means = matrix(c(h2*offspring_epsilon_o, rep(0.5*h2*offspring_epsilon_o,s + 2)), ncol = 3 + s )
  offspring_genetic_sibling_parents_eo = rmvnorm(n = sample_size, sigma = cov_matrix_parents_eg_sibs) + offspring_genetic_sibling_parents_eo_means
  if (s > 0) {
    colnames(offspring_genetic_sibling_parents_eo) = c("child_raw", sib.col.names, "dad", "mom")
  }
  #second set of sims are done here. No need to simulate environment, since it is implicitly done through the full lia.
  
  if ( s > 0) {
    T_val_sib = matrix(NA, ncol = 2, nrow = s)
    for (ii in 1:s) {
      cur_sib = paste("sib", ii,"_status", sep = "")
      T_val_sib[ii,] = thres[[1]][[status.cols[cur_sib][[1]] + 1]]
    }
    
    sibs.stat = unlist(status.cols[grep("sib", names(status.cols))])
    fam_cond_sibs_status = cond_sibs(sibs.stat = sibs.stat, h2 = h2, T_val_sib = T_val_sib, sample_size = sample_size)
    no.rows = 2*sample_size + dim(fam_cond_sibs_status)[1]
  } else {
    no.rows = 2*sample_size
  }
  
  
  #the holder of stuff:
  df = data.frame("child_raw" = numeric(no.rows))
  if (s > 0) {
    #create the columns needed for the siblings:
    for (i in 1:s) {
      df[[paste("sib", i, sep = "")]] = numeric(no.rows)
    }
    
  }
  df[["dad"]] = numeric(no.rows) 
  df[["mom"]] = numeric(no.rows)
  df[["child"]] = numeric(no.rows) 
  
  #### Storing simulated values ####
  
  #note the order is important, since we want the observations to match each other.
  #e.g. the sims with offspring being a case needs to match
  
  #genetic liability for offspring:
  df$child_raw = c(offspring_genetic_sibling_eps[,1],
                   offspring_genetic_sibling_parents_eo[,1],
                   fam_cond_sibs_status$child_raw)
  
  fam_cond_sibs =  unname(fam_cond_sibs_status[,grep("sib", colnames(fam_cond_sibs_status))])
  names(fam_cond_sibs) = NULL
  #Liability for siblings:
  df[,paste("sib", 1:s, sep = "")] = rbind(offspring_genetic_sibling_eps[,2:(s + 1)],
                                           offspring_genetic_sibling_parents_eo[,2:(s + 1)],
                                           fam_cond_sibs_status[,grep("sib", colnames(fam_cond_sibs_status))])
  
  #Full liability for offspring:
  df$child = c(offspring_genetic_sibling_eps[,1] + rnorm(n = sample_size, mean = 0, sd = sqrt(1 - h2)),
               offspring_epsilon_o,
               fam_cond_sibs_status$child_raw + rnorm(n = sample_size, mean = 0, sd = sqrt(1 - h2)))
  
  #Full liability for dad:
  df$dad = c(dad_cond,
             offspring_genetic_sibling_parents_eo[,(s + 2)],
             fam_cond_sibs_status$dad)
  #Full liability for mom:
  df$mom = c(mom_cond,
             offspring_genetic_sibling_parents_eo[,(s + 3)],
             fam_cond_sibs_status$mom)
  
  
  df = assign_grps_resamp(input = df, thresholds = thres, s = s, no.aoo.grps = no.aoo.grps)
  
  # cols = grep(paste(c( "sex", "status"), collapse = "|"), colnames(df))
  grp.means.resamp = df %>%
    group_by_at(vars(cols)) %>%
    summarise("obs" = list(child_raw)) %>%
    ungroup()
  return(grp.means.resamp)
}

cond_sibs = function(sibs.stat, h2, T_val_sib, sample_size) {
  s = length(sibs.stat)
  rho = 0.5*h2
  
  if (s == 1) {
    sib = rtruncnorm(n = sample_size, a = T_val_sib[1,1], b = T_val_sib[1,2])
    cov_cond_sibs = (diag(1 - (0.5*h2)^2 - 0.5*h2*(1 - 0.5*h2),3,3)) + 0.5*h2*(1 - 0.5*h2)
    cov_cond_sibs[1,1] = h2*(1 - h2)
    cov_cond_sibs[2, 3] = cov_cond_sibs[3,2]  = -1*(0.5*h2)^2
    fam_cond_sib = rmvnorm(n = sample_size, mean = rep(0,3), sigma = cov_cond_sibs) + matrix(rho*rep(sib, each = 3), byrow = T, ncol = 3)
    no.rows = sample_size
    df_fam = data.frame("child_raw" = numeric(no.rows),
                        "dad" = numeric(no.rows),
                        "mom" = numeric(no.rows),
                        "sib1" = numeric(no.rows))
    df_fam[,c("child_raw", "dad", "mom")] = fam_cond_sib
    df_fam[,"sib1"] = sib
    return(df_fam)
  }
  if (s == 2) {
    sib1 = rtruncnorm(n = sample_size, a = T_val_sib[1,1], b = T_val_sib[1,2])
    sib2 = rtruncnorm(n = sample_size, a = T_val_sib[2,1], b = T_val_sib[2,2])
    
    cov_cond_sibs = (diag(1 - (0.5*h2)^2 - 0.5*h2*(1 - 0.5*h2),3,3)) + 0.5*h2*(1 - 0.5*h2)
    cov_cond_sibs[1, 1] = h2*(1 - h2)
    cov_cond_sibs[2, 3] = cov_cond_sibs[3,2]  = -1*(0.5*h2)^2
    
    fam_cond_sib = rmvnorm(n = sample_size, mean = rep(0,3), sigma = cov_cond_sibs) + matrix(rho*rep((sib1 + sib2), each = 3), byrow = T, ncol = 3)
    no.rows = sample_size
    df_fam = data.frame("child_raw" = numeric(no.rows),
                        "dad" = numeric(no.rows),
                        "mom" = numeric(no.rows),
                        "sib1" = numeric(no.rows),
                        "sib2" = numeric(no.rows))
    df_fam[,c("child_raw", "dad", "mom")] = fam_cond_sib
    df_fam[,"sib1"] = sib1
    df_fam[,"sib2"] = sib2
    return(df_fam)
  }
  if (s >= 3) {
    cov_cond_sibs = (diag(1 - (0.5*h2)^2 - 0.5*h2*(1 - 0.5*h2),3,3)) + 0.5*h2*(1 - 0.5*h2)
    cov_cond_sibs[1, 1] = h2*(1 - h2)
    cov_cond_sibs[2, 3] = cov_cond_sibs[3,2]  = -1*(0.5*h2)^2
    trunc_sibs = list()
    for (ii in 1:s) {
      trunc_sibs[[paste("sib", ii, sep = "")]] = rtruncnorm(n = sample_size, a = T_val_sib[ii,1], b = T_val_sib[ii,2])
    }
    trunc_sibs = as.data.frame(trunc_sibs)
    mu_trunc_sibs = rowSums(trunc_sibs) * rho * (1 - rho) / (1 + rho - 2*rho^2)
    fam_cond_sibs = rmvnorm(n = sample_size, mean = rep(0,3), sigma = cov_cond_sibs) + matrix(rep(mu_trunc_sibs, each = 3), byrow = T, ncol = 3)
    no.rows = sample_size
    df_fam = data.frame("child_raw" = numeric(no.rows),
                        "dad" = numeric(no.rows),
                        "mom" = numeric(no.rows))
    df_fam[,c("child_raw", "dad", "mom")] = fam_cond_sibs
    for (ii in 1:s) {
      df_fam[[paste("sib",ii, sep = "")]] = trunc_sibs[[ii]]
    }
    return(df_fam)
  }
}

create.checker.mat = function(thres, s, no.age.grps = NULL){
  #"offspring_age" = 1:no.age.grps,
  expand.list = list("child_status" =  1:no.age.grps - 1)
  if (s > 0) {
    for (i in 1:s) {
      expand.list[[paste("sib", i, "_status", sep = "")]] = 1:no.age.grps - 1
    }
  }
  #  expand.list[["parent_age"]] = 2:no.age.grps # start set to 2, since we dont want parents and offsping with the same age(except for the last grp.)
  expand.list[["dad_status"]] = 1:no.age.grps - 1
  expand.list[["mom_status"]] = 1:no.age.grps - 1
  res = expand.grid(rev(expand.list))
  return(rev(res)) #rev used here to match the ordering of the grps with the output of dplyr
}

thres.assign.resamp = function(input, thresholds, s, no.age.grps, new.col.name, lia.col.name) {
  input[[new.col.name]] = 0
  for (jj in 1:no.age.grps) {
    input[[new.col.name]][input[[lia.col.name]] %between% thresholds[[1]][[jj]]] = jj - 1  
  }
  return(input)
}

assign_grps_resamp = function(input, thresholds, s, no.aoo.grps = NULL) {
  ## Need to simulate age and aoo_grp for the simulated data, and i need to account for it ..... 
  #do something else ?
  
  input = thres.assign.resamp(input = input, thresholds = thresholds, s = s, no.age.grps = no.age.grps,
                              new.col.name = "child_status", lia.col.name = "child")
  
  if (s > 0) {
    for (i in 1:s) {
      cur.sib = paste("sib", i, sep = "")
      cur.sib.stat = paste(cur.sib, "_status", sep = "")
      #     cur.sib.age = paste(cur.sib, "_age", sep = "")
      input = thres.assign.resamp(input = input, thresholds = thresholds, s = s, no.age.grps = no.age.grps,
                                  new.col.name = cur.sib.stat, lia.col.name = cur.sib)
    }
  }
  #currently using offsprings age + 1, and capping it at the highest age grp.
  input = thres.assign.resamp(input = input, thresholds = thresholds, s = s, no.age.grps = no.age.grps,
                              new.col.name = "dad_status", lia.col.name = "dad")
  input = thres.assign.resamp(input = input, thresholds = thresholds, s = s, no.age.grps = no.age.grps,
                              new.col.name = "mom_status", lia.col.name = "mom")
  return(input)
}

thres.assign = function(input, thresholds, s, no.age.grps, new.col.name, lia.col.name) {
  input[[new.col.name]] = 0
  for (ii in 1:no.age.grps) {
    input[[new.col.name]][input[[lia.col.name]] %between% thresholds[[1]][[ii]]] = ii - 1  
  }
  return(input)
}

assign_grps = function(input, thresholds, s, no.age.grps = NULL) {
  ## Need to simulate age and aoo_grp for the simulated data, and i need to account for it ..... 
  #do something else ?
  
  input = thres.assign(input = input, thresholds = thresholds, s = s, no.age.grps = no.age.grps, 
                       new.col.name = "child_status", lia.col.name = "child")
  
  if (s > 0) {
    for (i in 1:s) {
      cur.sib = paste("sib", i, sep = "")
      cur.sib.stat = paste(cur.sib, "_status", sep = "")
      #     cur.sib.age = paste(cur.sib, "_age", sep = "")
      input = thres.assign(input = input, thresholds = thresholds, s = s, no.age.grps = no.age.grps,
                           new.col.name = cur.sib.stat, lia.col.name = cur.sib)
    }
  }
  input = thres.assign(input = input, thresholds = thresholds, s = s, no.age.grps = no.age.grps,
                       new.col.name = "dad_status", lia.col.name = "dad")
  input = thres.assign(input = input, thresholds = thresholds, s = s, no.age.grps = no.age.grps,
                       new.col.name = "mom_status", lia.col.name = "mom")
  return(input)
}

sampler = function(thres, h2, s){
  no.age.grps = length(thres)
  cov_matrix = diag(1 - 0.5*h2, s + 3, s + 3) + 0.5*h2 
  cov_matrix[1,1] = h2
  cov_matrix[(s + 2),(s + 3)] = cov_matrix[(s + 3),(s + 2)] = 0
  
  child.env = rnorm(n = 500000, sd = sqrt(1 - h2))
  sim.dat = as.data.frame(rmvnorm(n = 500000, mean = rep(0, 3 + s), sigma = cov_matrix))
  if (s > 0) {
    colnames(sim.dat) = c("child_raw", paste("sib", 1:s, sep = ""), "dad", "mom")
  } else {
    colnames(sim.dat) = c("child_raw", "dad", "mom")
  }
  #filling out values for the "child" column, i.e. full child liability.
  sim.dat[["child"]] = sim.dat$child_raw + child.env 
  #move columns to be in the same order as the grp means columns:
  sim.dat[,2:(4 + s)] = sim.dat[,c("child", colnames(sim.dat)[-grep("child", colnames(sim.dat))])]
  #moving names to reflect the columns:
  colnames(sim.dat)[2:(4 + s)] = c("child", colnames(sim.dat)[-grep("child", colnames(sim.dat))])
  
  #sim.dat = assign_grps(input = sim.dat, thresholds = thres, s = s, no.age.grps = no.age.grps)
  
  #cols = colnames(sim.dat)[grep(paste(c("status", "age"), collapse = "|"), colnames(sim.dat))]
  #cur.means = sim.dat %>%
  #  group_by_at(vars(cols)) %>%
  #  summarise("obs" = list(child_raw)) %>%
  #  ungroup()
  #cur.means.str = apply(cur.means[,1:(3 + s)], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  grp.means = create.checker.mat(thres = thres, s, no.age.grps = no.age.grps)
  
  grp.thres = apply(as.matrix(grp.means[,1:(3 + s)]), 1, FUN = function(x) thres[x + 1]) # reversed dimensions
  
  grp.means = as_tibble(grp.means)
  #grp.means.str = apply(grp.means, MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  
  grp.means$obs <- list(NULL)
  checker.dat = sim.dat[,2:(4 + s)]
  for (i in 1:dim(grp.means)[1]) {
    grp.means$obs[[i]] = sim.dat$child_raw[rowSums(checker.dat < grp.thres[,i]) == 3 + s] #3 + s is os, siblings, and parents.
  }
  #  fill.entries = match(cur.means.str, grp.means.str)
  #  for (i in seq_along(fill.entries)) {
  #    ii = fill.entries[i]
  #    grp.means$obs[[ii]] = c(grp.means$obs[[ii]], cur.means$obs[[i]])
  #  }
  
  
  sem.vec = sapply(grp.means$obs, FUN = sem)
  
  
  
  nthreads = 10
  cl = makeCluster(nthreads, type = "SOCK")
  registerDoSNOW(cl)
  
  while (any(sem.vec > 0.01)) {
    cat("configurations requiring resampling:" , sum(sem.vec > 0.01), "\n"  )
    # ind.max.sem = which.max(sem.vec) ## maybe only take one at a time and update inbetween instead of all at once ?
    entries = sample(which(sem.vec > 0.01), size = nthreads)
    # grp.means[entries,]
    # age.cols = grp.means[ind.max.sem, str_subset(colnames(grp.means), "age")]
    status.cols = grp.means[entries, grep("status", colnames(grp.means))]
    
    ph = foreach(i = 1:nthreads,
                 .packages = c("data.table", "truncnorm", "mvtnorm", "dplyr")) %dopar% {
                   resampler(thres = thres, s = s, h2 = h2, status.cols = status.cols[i,], cols = cols)
                 }
    for (ii in seq_along(ph)) {
      update_entries = apply(ph[[ii]][,1:(3 + s)], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
      #    match.ph = match.ph[!is.na(match.ph)]
      if (!all(grp.means.str[grp.means.str %in% update_entries] == update_entries)) {
        print("STOOOOOP")
      }
      update_entries_no = match(update_entries,grp.means.str)
      grp.means$obs[update_entries_no] = map2(grp.means$obs[update_entries_no], ph[[ii]]$obs, ~ c(.x, .y))
      
    }
    
    #    grp.means[["mean"]] = map2_dbl(grp.means$all_obs, grp.means$all_resamp, ~ mean(c(.x,.y)))
    sem.vec = sapply(grp.means$obs, FUN = sem)
  }
  stopCluster(cl)
  grp.means[["means"]] = sapply(grp.means$obs, FUN = mean)
  
  return(grp.means)
}






