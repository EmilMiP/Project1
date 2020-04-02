setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("aoo_sampler.R")

assign_pheno = function(input,
                        h2 = 0.5,
                        child_aoo_col = "status",
                        age_col = "age",
                        dad_aoo_col = "dad_status",
                        mom_aoo_col = "mom_status",
                        age_parent = "age_parents",
                        sib_aoo_col_prefix = "sib",
                        sib_aoo_col_suffix = "_status",
                        new_col_name = "gen_lia",
                        nthreads = 10,
                        prev = .05,
                        no.grps = 2,
                        no.age.grps = 3,
                      #  thres,
                        s){
  
  ##### A much better way of dealing with the different grps is needed. we want the continuousness of the age variable combined with the discreteness of the configurations. some conversion to discrete and back again with interpolation ?
  prev.thres = prev*no.age.grps:1/no.age.grps
  thres = c(qnorm(1 - prev.thres))
  extreme_vals = c(round(min(input[,grep("age", colnames(input))])), round(max(input[,grep("age", colnames(input))])))
  age_cutoffs = quantile(seq(extreme_vals[1], extreme_vals[2]), prob = 1:(no.age.grps - 1)/no.age.grps)
  

  if (s > 0) {
    indivs = c("child", paste("sib", 1:s, sep = ""), "dad", "mom")
  } else {
    indivs = c("child", "dad", "mom")
  }
  
  for (i in seq_along(indivs)) {
    if (no.age.grps > 1) {
      input[[paste(indivs[i], "_status", sep = "")]] =  as.numeric(input[[paste(indivs[i], "_aoo", sep = "")]] > 0) + no.grps * rowSums(outer(input[[paste(indivs[i], "_age", sep = "")]], age_cutoffs, ">"))
    } else {
      input[[paste(indivs[i], "_status", sep = "")]] =  as.numeric(input[[paste(indivs[i], "_aoo", sep = "")]] > 0)
    }
  }
  input_str = apply(input[,grep("status", colnames(input))], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  actual.combs = unique(input_str)
  
  #assigning phenotype to people with aoo different than 0:
  input[[new_col_name]] = NA
  #estimating the genetic lia for controls:
  gen_lia_est = sampler(thres = thres, h2 = h2, s = s, actual.combs)

  gen_str = apply(gen_lia_est[,grep("grp", colnames(gen_lia_est))], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
#  data[[new_col_name]] = NA
  for (config in 1:length(gen_str)) {
    input[[new_col_name]][input_str %in% gen_str[config]] = gen_lia_est[config, "means"][[1]]
  }
  aoo_col = paste(indivs[1], "_aoo", sep = "")
  aoo.people = input[[aoo_col]] != 0
  #adding an artificial observation so the emperical cdf does not return 1 for the largest obs, meaning no Inf obs are generated.
  cum_prev = ecdf(c(input[[aoo_col]][aoo.people], max(input[[aoo_col]][aoo.people]) + 0.01 )) 
  #assigning the phenotype for cases, cased on their age of onset:
  input[[new_col_name]][aoo.people] = qnorm(1 - prev * (1 - cum_prev(input[[aoo_col]][aoo.people])), sd = sqrt(0.5*h2)) #data[[child_aoo_col]] -> input[["aoo"]]
  
  return(input)
}

options(scipen = 999)
library(data.table)
library(stringr)

#get list of all the summary stats to use, and their true values.
fileRoot = "/home/emp/faststorage/project1/simulatedData/"
#fileRoot = "D:\\Work\\Project1\\SimlatedData\\"
base = "\\\\AOOsibs2"
sword = "10k2"
cword = "C200"
bword = "NF"
identifier = ""
#exword = 


phen.files = list.files(path = fileRoot, pattern = "\\.phen", full.names = T)
phen.files = str_subset(phen.files, pattern = base)
phen.files = str_subset(phen.files, pattern = sword)
phen.files = str_subset(phen.files, pattern = cword)
phen.files = str_subset(phen.files, pattern = bword)
phen.files = str_subset(phen.files, pattern = identifier)


for (i in seq_along(phen.files)) {
  print(i)
  phen = as.data.frame(fread(phen.files[i]))
  ph = assign_pheno(input = phen, s = 2, no.age.grps = 1, no.grps = 2)
  fwrite(ph, phen.files[i], sep = " ", quote = F)
}


phen = as.data.frame(fread(paste(fileRoot, "AOOsibs2_10k2_NF_C200_v1.phen", sep = "")))
true = as.data.frame(fread(paste(fileRoot, "sibs2_10k2_NF_C200_v1.true", sep = "")))
ltfh = as.data.frame(fread(paste(fileRoot, "LTFHsibs2_10k2_NF_C200_v1.phen", sep = "")))
removers = as.data.frame(fread(paste(fileRoot, "sibs2_10k2_NF_C200_v1.excludeList", sep = "")))


ph = assign_pheno(input = phen, s = 2)
par(mfrow = c(2,2))
cor(true$offspringgeno_lia, ph$gen_lia)
plot(true$offspringgeno_lia, ph$gen_lia)
abline(b = 1, a = 0, col = "blue")

cor(true$offspringgeno_lia, ltfh$ltfh)
plot(true$offspringgeno_lia, ltfh$ltfh)
abline(b = 1, a = 0, col = "blue")


##
cor(true$offspringgeno_lia[-removers$FID], ph$gen_lia[-removers$FID])
plot(true$offspringgeno_lia[-removers$FID], ph$gen_lia[-removers$FID])
abline(b = 1, a = 0, col = "blue")

cor(true$offspringgeno_lia[-removers$FID], ltfh$ltfh[-removers$FID])
plot(true$offspringgeno_lia[-removers$FID], ltfh$ltfh[-removers$FID])
abline(b = 1, a = 0, col = "blue")

par(mfrow = c(1,1))



plot(ltfh$ltfh,ph$gen_lia)
abline(b = 1, a = 0, col = "blue")


c(0.4471322, 0.4482737, 0.7199584, 0.7159489 ) #3 grps
c(0.4462245, 0.4482737, 0.7196258, 0.7159489 ) # two grps
c(0.4452529, 0.4482737, 0.7191067, 0.7159489 ) # onegrp
c(0.4500616, 0.4482737, 0.7209058, 0.7159489 ) #one grp, halfheri
