library(truncnorm)
library(mvtnorm)
library(dplyr)
library(purrr)
library(stringr)
library(data.table)
library(doSNOW)
library(foreach)
library(parallel)
library(progress)

sem <- function(x){
  if (length(x) > 1) {
    return(sd(x)/sqrt(length(x)))
  } else {
    return(Inf)
  }
}

resampler = function(sample_size = 5e5, thres = thres, s = s, h2 = h2, status.cols, cols = cols, no.age.grps) {
  T_val_dad = thres[status.cols$dad_status + 1,]
  T_val_mom = thres[status.cols$mom_status + 1,]
  T_val_child = thres[status.cols$child_status + 1,]
  
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
      T_val_sib[ii,] = thres[status.cols[[cur_sib]] + 1,]
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
  
  
  df = assign_grps_resamp(input = df, thresholds = thres, s = s, no.age.grps = no.age.grps)
  
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

create.checker.mat = function(thres, s, no.age.grps = NULL, no.status = NULL){
  tot.grps = no.age.grps * no.status

  expand.list = list("child_grp" = 1:(tot.grps) - 1)
  if (s > 0) {
    for (i in 1:s) {
      expand.list[[paste("sib", i, "_grp", sep = "")]] = 1:(tot.grps) - 1
    }
  }
  #  expand.list[["parent_age"]] = 2:no.age.grps # start set to 2, since we dont want parents and offsping with the same age(except for the last grp.)
  expand.list[["dad_grp"]] = 1:(tot.grps) - 1
  expand.list[["mom_grp"]] = 1:(tot.grps) - 1
  
  res = expand.grid(rev(expand.list))
  return(rev(res)) #rev used here to match the ordering of the grps with the output of dplyr
}

thres.assign.resamp = function(input, thresholds, s, no.age.grps, new.col.name, lia.col.name) {
  input[[new.col.name]] = 0
  for (jj in 1:no.age.grps) {
    input[[new.col.name]][input[[lia.col.name]] %between% thresholds[jj,]] = jj - 1  
  }
  return(input)
}

assign_grps_resamp = function(input, thresholds, s, no.age.grps = NULL) {
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
#    input[[new.col.name]][input[[lia.col.name]] %between% thresholds[[1]][[ii]]] = ii - 1  
    input[[new.col.name]][input[[lia.col.name]] %between% thresholds[ii,]] = ii - 1  
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

sampler = function(thres, h2, s, actual.combs){
  cat("Drawing Initial samples. \n")
  no.age.grps = length(thres)
  cov_matrix = diag(1 - 0.5*h2, s + 3, s + 3) + 0.5*h2 
  cov_matrix[1,1] = h2
  cov_matrix[(s + 2),(s + 3)] = cov_matrix[(s + 3),(s + 2)] = 0
  
  child.env = rnorm(n = 5000000, sd = sqrt(1 - h2))
  sim.dat = as.data.frame(rmvnorm(n = 5000000, mean = rep(0, 3 + s), sigma = cov_matrix))
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
  
  cat("Assigning each individual to a group: \n")
  seq_n <- seq_len(nrow(sim.dat))
  all_split <- lapply(2:(3 + 1 + s), function(k) {
    print(k)
    if (length(thres) > 1) {
      sapply(1:no.age.grps, function(kk) {

        grp <- rowSums(outer(sim.dat[, k], thres[kk], `>`)) + no.age.grps * (kk - 1)
        split(seq_n, grp)
      })
    } else {
      grp <- rowSums(outer(sim.dat[, k], thres, `>`))
      split(seq_n, grp)
    }
  })
  
  cat("Assigning simulated family to a group: \n")
  checker.mat = as.matrix(create.checker.mat(thres, s = 2, no.age.grps = length(thres), no.status = 2))
  checker.mat = checker.mat[apply(checker.mat, MARGIN = 1, FUN = function(x) paste(x, collapse = "")) %in% actual.combs,]
  cat("")
  grp.cols = grep("grp" ,colnames(checker.mat))
  grp.means <- tibble::as_tibble(checker.mat)
  grp.means$obs <- lapply(1:nrow(grp.means), function(i) {
    tmp <- checker.mat[i, grp.cols]
    all_ind <- lapply(seq_along(tmp), function(k) {
      all_split[[k]][[tmp[[k]] + 1]]
    })
    sim.dat$child_raw[Reduce(intersect, all_ind)]
  })
  grp.means[["means"]] = sapply(grp.means$obs, FUN = mean)
  
  return(grp.means)
}





