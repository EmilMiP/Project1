library(data.table)
library(stringr)
library(ggplot2)

#get list of all the summary stats to use, and their true values.
#fileRoot = "/home/emp/faststorage/project1/simulatedData"
fileRoot = "D:\\Work\\Project1\\SimulatedData"
base = "sibs2"
sword = "100kx10k"
cword = "C1000"
bword = "_NF"
identifier = "maf005_prev0"
#exword = 
trueIdentifier = "maf005"


true.files = list.files(path = fileRoot, pattern = "\\.true", full.names = T)
true.files = str_subset(true.files, pattern = base)
true.files = str_subset(true.files, pattern = sword)
true.files = str_subset(true.files, pattern = cword)
true.files = str_subset(true.files, pattern = bword)
true.files = str_subset(true.files, pattern = trueIdentifier)
#true.files = str_subset(true.files, pattern = "ds")
#true.files = true.files[grep(true.files, pattern = cword)]
#true.files = true.files[-grep(true.files, pattern = exword)]
#true.files = true.files[grep(true.files, pattern = "10male2\\.t")]


#path to directory to find files in:
#bedRoot = "/home/emp/faststorage/project1/results/GWAS"
bedRoot = "D:\\Work\\Project1\\SimulatedData"
#identifying word for files to look for.

gwas.files = list.files(path = bedRoot, pattern = "CHILD_STATUS", full.names = T)
gwas.files = str_subset(gwas.files, base)
gwas.files = str_subset(gwas.files, bword)
gwas.files = str_subset(gwas.files, cword)
gwas.files = str_subset(gwas.files, sword)
gwas.files = str_subset(gwas.files, identifier)
gwas.files = gwas.files[-grep("\\.log$", gwas.files)] #remove log files
#gwas.files = gwas.files[grep(cword, gwas.files)]
#gwas.files = gwas.files[-grep(exword, gwas.files)]
#gwas.files = gwas.files[-grep("10male2_rev", gwas.files)]

ltfh.files = list.files(path = bedRoot, pattern = "_ltfh", full.names = T)
ltfh.files = str_subset(ltfh.files, base)
ltfh.files = str_subset(ltfh.files, bword)
ltfh.files = str_subset(ltfh.files, cword)
ltfh.files = str_subset(ltfh.files, sword)
ltfh.files = str_subset(ltfh.files, identifier)
ltfh.files = ltfh.files[-grep("\\.log$", ltfh.files)]
#ltfh.files = ltfh.files[-grep("new_ltfh", ltfh.files)]
#ltfh.files = ltfh.files[grep(cword, ltfh.files)]
#ltfh.files = ltfh.files[-grep(exword, ltfh.files)]
#ltfh.files = ltfh.files[-grep("10male2_rev", ltfh.files)]

ltfhNew.files = list.files(path = bedRoot, pattern = "omni_", full.names = T)
ltfhNew.files = str_subset(ltfhNew.files, base)
ltfhNew.files = str_subset(ltfhNew.files, bword)
ltfhNew.files = str_subset(ltfhNew.files, cword)
ltfhNew.files = str_subset(ltfhNew.files, sword)
#ltfhNew.files = str_subset(ltfhNew.files, identifier)
ltfhNew.files = ltfhNew.files[-grep("\\.log$", ltfhNew.files)]
ltfhNew.files = str_subset(ltfhNew.files, "/sibs2")


#ltfhNew.files = ltfhNew.files[grep(cword, ltfhNew.files)]
#ltfhNew.files = ltfhNew.files[-grep(exword, ltfhNew.files)]
#ltfhNew.files = ltfhNew.files[-grep("10male2_rev", ltfhNew.files)]

gwax.files = list.files(path = bedRoot, pattern = "_GWAX", full.names = T)
gwax.files = str_subset(gwax.files, base)
gwax.files = str_subset(gwax.files, bword)
gwax.files = str_subset(gwax.files, cword)
gwax.files = str_subset(gwax.files, sword)
gwax.files = str_subset(gwax.files, identifier)
gwax.files = gwax.files[-grep("\\.log$", gwax.files)]

#gwax.files = gwax.files[grep(cword, gwax.files)]
#gwax.files = gwax.files[-grep(exword, gwax.files)]
#gwax.files = gwax.files[-grep("10male2_rev", gwax.files)]

true.gwas.files = list.files(path = bedRoot, pattern = "offspringgeno", full.names = T)
true.gwas.files = str_subset(true.gwas.files, base)
true.gwas.files = str_subset(true.gwas.files, bword)
true.gwas.files = str_subset(true.gwas.files, cword)
true.gwas.files = str_subset(true.gwas.files, sword)
true.gwas.files = str_subset(true.gwas.files, identifier)
true.gwas.files = true.gwas.files[-grep("\\.log$", true.gwas.files)] #remove log files
#true.gwas.files = str_subset(true.gwas.files, "NF")
#true.gwas.files = true.gwas.files[-grep("10male2_rev", true.gwas.files)]

all.files = c(gwas.files,gwax.files, ltfh.files, ltfhNew.files,true.gwas.files)
version.subset = list()
pattern.vec = paste("v", c(1,10,2:9), "_", sep = "")
for (i in 1:length(true.files)) {
  version.subset[[i]] = str_subset(all.files, pattern.vec[i])
}

version.len = length(version.subset[[1]])
version.names = sub( ".*v1_", "", version.subset[[1]])
#version.names = gsub("CHILD_STATUS", "GWAS", version.names)
#version.names = gsub("gendered_stat", "Newltfh_nosib", version.names)
#version.names = gsub("GENDERED_", "", version.names)
#version.names = sort(version.names)

#true vals
val.list = list()
for (i in 1:length(true.files)) {
  #true values:
  datrue = fread(true.files[i])
  lia_beta = datrue$lia_betas[1:10000] #currently hardcoded for M = 10k.
  colnames(datrue) = gsub("\\.", "_", colnames(datrue))
  mat = matrix(NA, ncol = 4, nrow = version.len)
  colnames(mat) = c("no. false pos.","no Sig Snps", "chisq_null", "chisq_causal")
  rownames(mat) = version.names
  
  for (j in 1:version.len) {
    #first pheno:
    da = fread(version.subset[[i]][j])
    da[["chisq"]] = qchisq(da$P, df = 1, lower.tail = FALSE)
    #proportion of true sigg snps
    mat[j,2] = sum((da$P <= 5e-08) & (lia_beta  != 0))
    mat[j,1] = sum((da$P <= 5e-08)) - mat[j,2]
    #mean of true non sigg snps
    mat[j,3] = mean(da$chisq[lia_beta == 0]) #we are just checking if they are causal
    #mean of true sigg snps
    mat[j,4] = mean(da$chisq[lia_beta != 0])
  }
  val.list[[i]] = mat
}

bs.res = matrix(NA, ncol = 8, nrow = version.len)
rownames(bs.res) = version.names
colnames(bs.res) = c("mean_false_pos.", "mean_false_pos._sd","mean_C_snps","mean_C_snps_sd","chisq_null_mean","chisq_null_sd", "chisq_causal_mean", "chisq_causal_sd")
for (i in 1:dim(bs.res)[1]) {
  data = matrix(unlist(lapply(val.list, FUN = function(x) x[i,])), ncol = 4, byrow = T) #extracts all obs for one method
  ctr = 1 #used to place obs and their standard deviations in the correct place
  for (j in 1:4) {
    bs.dat = replicate(10000, sample(data[,j], replace = T))
    bs.means = colMeans(bs.dat)
    bs.res[i,ctr] = mean(bs.means)
    bs.res[i,ctr + 1] = sd(bs.means)
    ctr = ctr + 2 #add two because we save two values per iteration.
  }
}
bs.res

bs.res = as.data.frame(bs.res)
bs.res[["prev"]] = c(rep(c(" 5% Prevalence", "10% Prevalence", "20% Prevalence"), 6), "true")

bs.res[["type"]] = c(rep(c("GWAS", "GWAX","LTFH"),each = 3), rep(c("NewLTFH_sex_1grp", "NewLTFH_sex_2grp", "NewLTFH_sex_5grp"), each = 3), "true")
bs.res[["order"]] = paste(bs.res$type, bs.res$prev, sep = "_")
#bs.res$order = factor(bs.res$order, levels = paste(bs.res$type, bs.res$offset, sep = "_"))
#bs.res = bs.res[order(bs.res$offset),]
#"ltfh_nosib",
#rep("ltfh_nosib", 6),
plot_dat = bs.res[-nrow(bs.res),]
c = datrue$No_Causal_SNPs[1]
p = ggplot(plot_dat, aes(x = type, y = mean_C_snps/c, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean_C_snps - mean_C_snps_sd) / c, ymax = (mean_C_snps + mean_C_snps_sd) / c)) +
  theme(axis.text.x = element_text(angle = -90)) +
  xlab("") +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(palette = "Blues") + 
  ylab("Fraction of causal SNPs detected") + 
  facet_wrap(~prev, ncol = 1)
plot(p)

chi_null_plot = ggplot(plot_dat, aes(color = type)) +
  geom_pointrange(aes(x = type, y = chisq_null_mean, ymin = chisq_null_mean - 1.96*(chisq_null_sd/sqrt(10)), ymax = chisq_null_mean + 1.96*(chisq_null_sd/sqrt(10)))) +
  geom_hline(yintercept = 1) +
#  geom_errorbar(aes(ymin = chisq_null_mean - 1.96*(chisq_null_sd/sqrt(10)) , ymax = chisq_null_mean + 1.96*(chisq_null_sd/sqrt(10)))) +
  xlab("") +
#  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(palette = "Blues") + 
  ylab("Mean Null Chisq") + 
#  ylim(min = 0.9, max = 1.1) +
  facet_wrap(~prev, ncol = 1)
plot(chi_null_plot)

chi_causal_plot = ggplot(plot_dat, aes(color = type)) +
  geom_pointrange(aes(x = type, y = chisq_causal_mean , ymin = chisq_causal_mean - 1.96*(chisq_causal_sd/sqrt(10)), ymax = chisq_causal_mean  + 1.96*(chisq_causal_sd/sqrt(10)))) +
  #  geom_errorbar(aes(ymin = chisq_null_mean - 1.96*(chisq_null_sd/sqrt(10)) , ymax = chisq_null_mean + 1.96*(chisq_null_sd/sqrt(10)))) +
  xlab("") +
  #  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(palette = "Blues") + 
  ylab("Mean Causal Chisq") +
  facet_wrap(~prev, ncol = 3)
plot(chi_causal_plot)






#loading and filtering
true = fread(true.files[1])
true.betas.pos = which(true$lia_betas[1:10000] != 0)
true.betas = true$lia_betas[true.betas.pos]
prevword   = "prev020"
prev.files = str_subset(version.subset[[1]], prevword)
#storage
beta_mat   = matrix(NA, ncol = length(prev.files), nrow = length(true.betas))
p_mat      = matrix(NA, ncol = length(prev.files), nrow = length(true.betas))
#names
temp = sapply(str_split(prev.files, "maf005_"), FUN = function(x) tail(x,n = 1))
colnames(beta_mat) = sapply(str_split(temp, pattern = "\\."), FUN = function(x) x[1])
colnames(p_mat)    = sapply(str_split(temp, pattern = "\\."), FUN = function(x) x[1])
#extraction
for (i in seq_along(prev.files)) {
  da = fread(prev.files[i])
  if (!any(colnames(da) == "BETA")) {
    beta_mat[,i] = log(da$OR[true.betas.pos])
  } else {
    beta_mat[,i] = da$BETA[true.betas.pos]
  }
  p_mat[,i] = -log10(da$P[true.betas.pos])
}
#constructing dataset for plotting
df = as.data.frame(matrix(NA, ncol = 3, nrow = nrow(beta_mat)*ncol(beta_mat)))
colnames(df) = c("beta.est", "method", "-log10(p)")
M = nrow(beta_mat)
for (i in 1:ncol(beta_mat)) {
  df[(1 + (i - 1) * M):(i * M),1] = beta_mat[,i]
  df[(1 + (i - 1) * M):(i * M),2] = rep(colnames(beta_mat)[i], M)
  df[(1 + (i - 1) * M):(i * M),3] = p_mat[,i]
}

df$true = true.betas
ggplot(df, aes(x = beta.est, y = true, alpha = .1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~method)

p_mat = as.data.frame(p_mat)

plotter = function(data, column1, column2) {
  p <- ggplot(data) +
       geom_point(aes(x = .data[[column1]] , y = .data[[column2]])) + 
       geom_hline( yintercept = -log10(5e-8)) +
       geom_vline( xintercept = -log10(5e-8))
  return(p)
}

plot_list = list()
ctr = 1
for (i in 2:ncol(p_mat)) {
  for (j in 1:(i-1)) {
    plot_list[[ctr]] = plotter(p_mat, colnames(p_mat)[i], colnames(p_mat)[j])
    ctr = ctr + 1
  }
}


library(egg)
ggarrange(plots = plot_list)

# pentioneret kode --------------------------------------------------------


#
#
#
#phen.files = list.files(path = fileRoot, pattern = "v1_", full.names = T)
#phen.files = phen.files[grep(cword, phen.files)]
#phen.files = phen.files[grep(base, phen.files)]
#phen.files = phen.files[-grep("map|\\.ped|exclude", phen.files)]
#
#true.dist = phen.files[grep("true", phen.files)]
#phen.dist = phen.files[-grep("true|LTFH", phen.files)]
#ltfh.dist = phen.files[grep("LTFH", phen.files)]
#
#
#true = fread(true.dist)
#phen = fread(phen.dist)
#phen[["lia"]] = true$offspring_lia
#phen[["true_gen_lia"]] = true$offspringgeno_lia
#
#for (i in seq_along(ltfh.dist)) {
#  identifier = substr(tail(strsplit(ltfh.dist[i], "\\\\")[[1]], n = 1), start = 5, stop = 8)
#  ltfh = fread(ltfh.dist[i])
#  phen[[paste("ltfh", identifier, sep = "")]] = ltfh$ltfh
#  if (identifier == "_5_5") {
#    phen = cbind(phen, ltfh[,-(1:5)])
#  }
#}
#phen = as.data.frame(phen)
#colnames(phen)[colnames(phen) == "CHILD_STATUS"] = "CHILD_STATUS_5_5"
#colnames(phen)[colnames(phen) == "GWAX"] = "GWAX_5_5"
#cols = colnames(phen)[grep("lia|CHILD_STATUS|ltfh|GWAX", colnames(phen))]
#cor(phen[,cols])[c("lia", "true_gen_lia"),]
#
#
#pdat = phen[,cols]
#pdat_5_5 = melt(pdat[,grep("_5_5|true", colnames(pdat))], id.vars = "true_gen_lia", variable.name = "method")
#p = ggplot(pdat_5_5, aes(x = true_gen_lia, y = value, color = method)) + 
#  geom_point() + 
#  geom_abline(intercept = 0, slope = 1) +
#  facet_grid(. ~ method) +
#  theme_minimal()
#plot(p)
#
#betas = true$lia_betas[1:10000] + true$male_betas[1:10000]
#
#
#files = str_subset(all.files, "v1_")
#files = str_subset(files, "_5_5")
#files = str_subset(files, "qassoc|logistic")
#sig.thres = 5e-8
#
#beta.names = sapply(str_split(files, pattern = "\\_10male_|\\."), FUN = function(x) x[2])
#beta.est = data.frame("true.betas" = betas)
#p.est = data.frame("true.betas" = betas)
#for (i in 1:length(files)) {
#  betas.ph = fread(files[i])
#  if ("BETA" %in% colnames(betas.ph)) {
#    beta.est[[beta.names[i]]] = betas.ph$BETA 
#    #  beta.est[[paste(beta.names[i],"_col", sep = "")]] = c("black","red")[(betas.ph$P <= sig.thres) + 1]
#  } else {
#    beta.est[[beta.names[i]]] = log(betas.ph$OR)
#    # beta.est[[paste(beta.names[i],"_col", sep = "")]] = c("black","red")[(betas.ph$P <= sig.thres) + 1]
#  }
#}
#beta.est = melt(beta.est, id.vars = "true.betas", measure.vars = grep(pattern = "_5_5",colnames(beta.est)),
#                variable.name = "method")
#beta.est[["sig"]] = (betas.ph$P <= sig.thres) + 1
#
#
#p = ggplot(beta.est, aes(x = true.betas, y = value,color = sig)) + 
#  geom_point() + 
#  geom_abline(intercept = 0, slope = 1) +
#  facet_grid(. ~ method) +
#  theme_minimal()
#plot(p)
#
#
#
#
#
##
#
#
#
#
#
#
#