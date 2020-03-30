options(scipen = 999)
library(data.table)
library(stringr)

#get list of all the summary stats to use, and their true values.
#fileRoot = "/home/emp/faststorage/project1/simulatedData"
fileRoot = "D:\\Work\\Project1\\SimlatedData\\"
base = "sibs2"
sword = "10k2"
cword = "C200"
bword = "F"
identifier = "_NDS"
#exword = 


true.files = list.files(path = fileRoot, pattern = "\\.true", full.names = T)
true.files = str_subset(true.files, pattern = base)
true.files = str_subset(true.files, pattern = sword)
true.files = str_subset(true.files, pattern = cword)
true.files = str_subset(true.files, pattern = bword)
true.files = str_subset(true.files, pattern = identifier)





k = .05
no.grps = 5
prev.thres = k*no.grps:1/no.grps
#prevs = c()
#for (i in 1:no.grps) {
#  prevs[i] = (prev.thres[i + 1] + prev.thres[i]) / 2
#}
#


#aoo.thres = as_tibble(matrix(list(), ncol = 1, nrow = no.grps ))
#colnames(aoo.thres) = paste("age_grp", 1:no.grps, sep = "")
#colnames(aoo.thres) = paste("grp", 0:no.grps, sep = 
aoo.thres = matrix(NA, ncol = 3, nrow = length(prev.thres))

for (j in 1:(length(prev.thres))) {
  #aoo.thres[[1]][j] = list(c(-Inf, qnorm(1 - prev.thres[j])))
  aoo.thres[j,] = c(-Inf, qnorm(1 - prev.thres[j]), Inf)
} 

aoo.thres[length(prev.thres),] = c(qnorm(1 - prev.thres[length(prev.thres)]), Inf)

thres = matrix(c(-Inf, rep(qnorm(1 - .05), 2), Inf), ncol = 2)





min.age = 30
span.of.years = 50


for (i in seq_along(true.files))  {
  print(i)
  true = fread(true.files[i])

  phen = data.frame("FID" = 1:dim(true)[1], "IID" = 1:dim(true)[1])
  
  cases = true$offspring_lia[true$offspring_lia >= qnorm(1 - .05)]
  #constructing an s-curved cumulative distribution function in the interval of cases:
  
# #data = rbeta(n = 10000, shape1 = 2, shape2 = 100)
# mindat = min(cases)
# maxdat = max(cases)
# norm.dat = (cases - mindat) / (maxdat - mindat)
# ph = min(cases) + norm.dat * (max(cases) - min(cases)) 
#  cum.prev = ecdf(c(min(ph) - .001,ph, max(ph) + .001)) #addeding a very minor increase to the max obs to prevent infinity later on.
  phen["child_aoo"] = 0
  #scaling aoo to be the same interval as the age:
  #transforming vals to be between 0 and 1, then scale up.
  phen$child_aoo[true$offspring_lia >= qnorm(1 - .05)] = (cases - min(cases)) / max(cases - min(cases)) * span.of.years + min.age   # cum.prev(cases) 
  new.subset = phen$child_aoo != 0
  phen["child_age"] = runif(n = dim(phen)[1], min = 0, max = span.of.years)
  s = 2
  for (ii in 1:s) {
    cur.sib = paste("sib", ii, "_aoo", sep = "")
    phen[[cur.sib]] = 0
    cur.sib.lia = paste("sibling_lia.", ii, sep = "")
    cases = true[[cur.sib.lia]][true[[cur.sib.lia]] >= qnorm(1 - .05)]
    phen[[cur.sib]][true[[cur.sib.lia]] >= qnorm(1 - .05)] = (cases - min(cases)) / max(cases - min(cases)) * span.of.years + min.age #  cum.prev(cases) #
    phen[[paste("sib", ii, "_age", sep = "")]] = phen$child_age
  }
  
  
  phen[["dad_aoo"]] = 0
  cases = true$parents_lia.1[true$parents_lia.1 >= qnorm(1 - .05)]
  phen$dad_aoo[true$parents_lia.1 >= qnorm(1 - .05)] = (cases - min(cases)) / max(cases - min(cases)) * span.of.years + min.age  #  cum.prev(cases) #
  phen[["dad_age"]] = phen$child_age + min.age

  phen[["mom_aoo"]] = 0
  cases = true$parents_lia.2[true$parents_lia.2 >= qnorm(1 - .05)]
  phen$mom_aoo[true$parents_lia.2 >= qnorm(1 - .05)] = (cases - min(cases)) / max(cases - min(cases)) * span.of.years + min.age  #  cum.prev(cases) #
  phen[["mom_age"]] = phen$child_age + min.age
  
  file.dist = paste(fileRoot,"AOO", base, "_", sword, "_", bword, "_", cword, "_v", i,identifier, ".phen", sep = "")
  fwrite(phen, file.dist, sep = " ", quote = FALSE)
}


phen.files = list.files(path = fileRoot, pattern = "\\.phen", full.names = T)
phen.files = str_subset(phen.files, pattern = base)
phen.files = str_subset(phen.files, pattern = sword)
phen.files = str_subset(phen.files, pattern = cword)
phen.files = str_subset(phen.files, pattern = bword)
phen.files = str_subset(phen.files, pattern = identifier)

#phen.files = str_subset(phen.files, pattern = "/sibs")

for (i in seq_along(phen.files)) {
  phen = fread(phen.files[i])
  phen[["GWAX"]] = as.numeric(rowSums(phen[,c("CHILD_STATUS", "P1_STATUS", "P2_STATUS", "SIB_STATUS")]) > 0)
  fwrite(phen, phen.files[i], sep = " ", quote = F)
}



#phen[["gen_lia"]] = NA
#frac.of.cases = cum.prev(true$offspring_lia[new.subset])
#h2 = .5
#phen$gen_lia[new.subset] = qnorm(1 - k * (1 - frac.of.cases), sd = sqrt(h2))
#phen$

#
#sq = function(x) {return(x^2)  }
#
#
#f = approxfun(1:10/2, sq(1:10/2))
#f(2.3)
#sq(2.3)
#
#phen3 = phen[new.subset,]
#ltfh2 = ltfh[new.subset,]
#cor(phen2$gen_lia[new.subset], true$offspringgeno_lia[new.subset])
#plot(true$offspringgeno_lia[new.subset], phen2$gen_lia[new.subset])
#abline(b = 1, a = 0)
#
#true2 = true$offspringgeno_lia[new.subset]
#par(mfrow = c(2,2))
#cor(true2, phen2$gen_lia[new.subset])
#plot(true2, phen2$gen_lia[new.subset])
#abline(b = 1, a = 0)
#
#cor(true2, phen3$aoo_ltfh)
#plot(true2, phen3$aoo_ltfh)
#abline(b = 1, a = 0)
#
#cor(true2, ltfh2$new_ltfh_5_5)
#plot(true2, ltfh2$new_ltfh_5_5)
#abline(b = 1, a = 0)
#
#cor(true2, ltfh2$ltfh)
#plot(true2, ltfh2$ltfh)
#abline(b = 1, a = 0)



ph = split(phen2, cut(phen2$age, min.age:(min.age + span.of.years)))
sapply(ph, FUN = function(x) sum(x$aoo != 0))
cumsum(sapply(ph, FUN = function(x) sum(x$aoo != 0))) / (sum(phen2$aoo != 0) +  1) # adding 1 here to avoid inf, difference occurs at 3rd-4th decimal





# parkeret ----------------------------------------------------------------



phen[["aoo"]] = 0
for (i in 1:no.grps) {
  phen$aoo[true$offspring_lia %between% c(qnorm(1 - prev.thres[i + 1]), qnorm(1 - prev.thres[i]))] = no.grps + 1 - i
  #print(sum(true$offspring_lia %between% c(qnorm(1 - prev.thres[i + 1]), qnorm(1 - prev.thres[i])) ) )
}

phen[["aoo_p1"]] = 0
phen[["aoo_p2"]] = 0
for (i in 1:no.grps) {
  phen$aoo_p1[true$parents_lia.1 %between% c(qnorm(1 - prev.thres[i + 1]), qnorm(1 - prev.thres[i]))] = no.grps + 1 - i
  phen$aoo_p2[true$parents_lia.2 %between% c(qnorm(1 - prev.thres[i + 1]), qnorm(1 - prev.thres[i]))] = no.grps + 1 - i
  #print(sum(true$offspring_lia %between% c(qnorm(1 - prev.thres[i + 1]), qnorm(1 - prev.thres[i])) ) )
}

phen[["aoo_sib1"]] = 0
phen[["aoo_sib2"]] = 0
for (i in 1:no.grps) {
  phen$aoo_sib1[true$sibling_lia.1 %between% c(qnorm(1 - prev.thres[i + 1]), qnorm(1 - prev.thres[i]))] = no.grps + 1 - i
  phen$aoo_sib2[true$sibling_lia.1 %between% c(qnorm(1 - prev.thres[i + 1]), qnorm(1 - prev.thres[i]))] = no.grps + 1 - i
  #print(sum(true$offspring_lia %between% c(qnorm(1 - prev.thres[i + 1]), qnorm(1 - prev.thres[i])) ) )
}

phen[["age"]] = sample(1:no.grps, size = dim(true)[1], replace = T)
phen[["age_parents"]] = phen$age + 1
phen$age_parents[phen$age_parents == 5] = 4

thres.assign = function(input, lia.dat, thresholds, s, no.age.grps, no.aoo.grps, new.col.name, lia.col.name, age.col.name) {
  input[[new.col.name]] = 0
  for (ii in 1:no.age.grps) {
    for (jj in 1:no.aoo.grps) {
      input[[new.col.name]][input[[age.col.name]] == ii & (lia.dat[[lia.col.name]] %between% thresholds[,ii][[1]][[jj]])] = jj - 1  
    }
  }
  return(input)
}

assign_aoo_grps = function(input, lia.dat, thresholds, s, no.age.grps = NULL, no.aoo.grps = NULL) {
  ## Need to simulate age and aoo_grp for the simulated data, and i need to account for it ..... 
  #do something else ?
  
  input = thres.assign(input = input, lia.dat = lia.dat, thresholds = thresholds, s = s, no.age.grps = no.age.grps, no.aoo.grps = no.aoo.grps,
                       new.col.name = "child_status", lia.col.name = "offspring_lia", age.col.name = "age")
  
  if (s > 0) {
    for (i in 1:s) {
      cur.sib = paste("sib", i, sep = "")
      cur.sib.stat = paste(cur.sib, "_status", sep = "")
      cur.sib.lia = paste("sibling_lia.", i, sep = "")
      #     cur.sib.age = paste(cur.sib, "_age", sep = "")
      input = thres.assign(input = input, lia.dat = lia.dat, thresholds = thresholds, s = s, no.age.grps = no.age.grps, no.aoo.grps = no.aoo.grps,
                           new.col.name = cur.sib.stat, lia.col.name = cur.sib.lia, age.col.name = "age")
    }
  }
  #currently using offsprings age + 1, and capping it at the highest age grp.
  input = thres.assign(input = input, lia.dat = lia.dat, thresholds = thresholds, s = s, no.age.grps = no.age.grps, no.aoo.grps = no.aoo.grps,
                       new.col.name = "dad_status", lia.col.name = "parents_lia.1", age.col.name = "age_parents")
  input = thres.assign(input = input, lia.dat = lia.dat, thresholds = thresholds, s = s, no.age.grps = no.age.grps, no.aoo.grps = no.aoo.grps,
                       new.col.name = "mom_status", lia.col.name = "parents_lia.2", age.col.name = "age_parents")
  return(input)
}




phen = assign_aoo_grps(input = phen, lia.dat = true, thresholds = aoo.thres, s = s, no.age.grps = 4, no.aoo.grps = 5)



table("age" = phen$age, "grp" = phen$aoo)
table("age" = phen$age, "grp" = phen$child_status)
phen %>% group_by(age, aoo_grp) %>% summarise(n())

#curve(expr = .5*(x^2 - .1^2), from = 0.1, to = 1)
phen = phen[,-grep(colnames(phen), pattern = "aoo")]
phen = phen[, c("FID", "IID", "age", "child_status", "sib1_status", "sib2_status", "age_parents", "dad_status", "mom_status")]
colnames(grp.means)
colnames(phen)


est_config = apply(grp.means[,1:(5 + s)], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
data_str = apply(phen[,1:(5 + s) + 2], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
phen[["aoo_ltfh"]] = NA
for (config in 1:length(est_config)) {
  phen[["aoo_ltfh"]][data_str %in% est_config[config]] = grp.means[config, "means"][[1]]
}

hist(phen$aoo_ltfh)
plot(true$offspringgeno_lia, phen$aoo_ltfh, pch = 19)
abline(b = 1, a = 0)
ltfh = fread("D:\\data\\workdata\\707293\\FIUN\\SimData\\pedFiles\\LTFH_5_5sibs2_10k2_F_C200_v1_NDS.phen")
excludelist = fread("D:\\data\\workdata\\707293\\FIUN\\SimData\\pedFiles\\sibs2_10k2_F_C200_v1_NDS.excludeList_5_5")
removers = !(phen$FID %in% excludelist$FID)
cor(true$offspringgeno_lia[removers,], cbind(phen$aoo_ltfh, ltfh$ltfh, ltfh$new_ltfh_5_5)[removers,])
