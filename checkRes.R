library(data.table)
library(stringr)
library(ggplot2)

#get list of all the summary stats to use, and their true values.
fileRoot = "D:\\data\\workdata\\707293\\FIUN\\SimData\\pedFiles\\"
base = "sibs2"
sword = "10k2"
cword = "C100_"
bword = "_NF"
identifier = "10male2_rev"
#exword = 



true.files = list.files(path = fileRoot, pattern = "\\.true", full.names = T)
true.files = str_subset(true.files, pattern = base)
true.files = str_subset(true.files, pattern = sword)
true.files = str_subset(true.files, pattern = cword)
true.files = str_subset(true.files, pattern = bword)
true.files = str_subset(true.files, pattern = identifier)
#true.files = str_subset(true.files, pattern = "ds")
#true.files = true.files[grep(true.files, pattern = cword)]
#true.files = true.files[-grep(true.files, pattern = exword)]
#true.files = true.files[grep(true.files, pattern = "10male2\\.t")]


#path to directory to find files in:
bedRoot = "D:\\data\\workdata\\707293\\FIUN\\SimData\\bedFiles\\"
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
ltfh.files = ltfh.files[-grep("new_ltfh", ltfh.files)]
#ltfh.files = ltfh.files[grep(cword, ltfh.files)]
#ltfh.files = ltfh.files[-grep(exword, ltfh.files)]
#ltfh.files = ltfh.files[-grep("10male2_rev", ltfh.files)]

ltfhNew.files = list.files(path = bedRoot, pattern = "new_ltfh", full.names = T)
ltfhNew.files = str_subset(ltfhNew.files, base)
ltfhNew.files = str_subset(ltfhNew.files, bword)
ltfhNew.files = str_subset(ltfhNew.files, cword)
ltfhNew.files = str_subset(ltfhNew.files, sword)
ltfhNew.files = str_subset(ltfhNew.files, identifier)
ltfhNew.files = ltfhNew.files[-grep("\\.log$", ltfhNew.files)]

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

true.gwas.files = list.files(path = bedRoot, pattern = "_lia", full.names = T)
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

version.len = length(version.subset[[3]])
version.names = sub( ".*v2_", "", version.subset[[3]])
version.names = gsub("CHILD_STATUS", "GWAS", version.names)
#version.names = gsub("gendered_stat", "Newltfh_nosib", version.names)
version.names = gsub("GENDERED_", "", version.names)
version.names = sort(version.names)

#true vals
val.list = list()
for (i in 1:length(true.files)) {
  #true values:
  datrue = fread(true.files[i])
  datrue = datrue[1:10000,]
  colnames(datrue) = gsub("\\.", "_", colnames(datrue))
  mat = matrix(NA, ncol = 4, nrow = version.len)
  colnames(mat) = c("no. false pos.","no Sig Snps", "chisq_null", "chisq_causal")
  rownames(mat) = version.names
  
  for (j in 1:version.len) {
    #first pheno:
    da = fread(version.subset[[i]][j])
    da[["chisq"]] = qchisq(da$P, df = 1, lower.tail = FALSE)
    #proportion of true sigg snps
    mat[j,2] = sum((da$P <= 5e-08) & (datrue$lia_betas_1  != 0))
    mat[j,1] = sum((da$P <= 5e-08)) - mat[j,2]
    #mean of true non sigg snps
    mat[j,3] = mean(da$chisq[datrue$lia_betas_1 == 0]) #we are just checking if they are causal
    #mean of true sigg snps
    mat[j,4] = mean(da$chisq[datrue$lia_betas_1 != 0])
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
#version.names = sub("GWAS", "GWAS_5_5", version.names)

bs.res = as.data.frame(bs.res)
bs.res[["offset"]] = c(rep(c("5_5", "1_9","2_8", "3_7", "4_6"), 2), rep(c("1_9","2_8", "3_7", "4_6", "5_5"), 2))
bs.res[["type"]] = rep(c("GWAS", "GWAX","LTFH","NewLTFH"),each = 5)
bs.res[["order"]] = paste(bs.res$type, bs.res$offset, sep = "_")
bs.res$order = factor(bs.res$order, levels = paste(bs.res$type, bs.res$offset, sep = "_"))
bs.res = bs.res[order(bs.res$offset),]
#"ltfh_nosib",
#rep("ltfh_nosib", 6),
c = datrue$No_Causal_SNPs[1]
p = ggplot(bs.res, aes(x = type, y = mean_C_snps/c, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean_C_snps - mean_C_snps_sd) / c, ymax = (mean_C_snps + mean_C_snps_sd) / c)) +
  theme(axis.text.x = element_text(angle = -90)) +
  xlab("") +
  coord_flip() +
  theme_minimal() +
  facet_wrap(~offset, ncol = 1)
plot(p)


phen.files = list.files(path = fileRoot, pattern = "v1_", full.names = T)
phen.files = phen.files[grep(cword, phen.files)]
phen.files = phen.files[grep(base, phen.files)]
phen.files = phen.files[-grep("map|\\.ped|exclude", phen.files)]

true.dist = phen.files[grep("true", phen.files)]
phen.dist = phen.files[-grep("true|LTFH", phen.files)]
ltfh.dist = phen.files[grep("LTFH", phen.files)]


true = fread(true.dist)
phen = fread(phen.dist)
phen[["lia"]] = true$offspring_lia
phen[["true_gen_lia"]] = true$offspringgeno_lia

for (i in seq_along(ltfh.dist)) {
  identifier = substr(tail(strsplit(ltfh.dist[i], "\\\\")[[1]], n = 1), start = 5, stop = 8)
  ltfh = fread(ltfh.dist[i])
  phen[[paste("ltfh", identifier, sep = "")]] = ltfh$ltfh
  if (identifier == "_5_5") {
    phen = cbind(phen, ltfh[,-(1:5)])
  }
}
phen = as.data.frame(phen)
colnames(phen)[colnames(phen) == "CHILD_STATUS"] = "CHILD_STATUS_5_5"
colnames(phen)[colnames(phen) == "GWAX"] = "GWAX_5_5"
cols = colnames(phen)[grep("lia|CHILD_STATUS|ltfh|GWAX", colnames(phen))]
cor(phen[,cols])[c("lia", "true_gen_lia"),]


pdat = phen[,cols]
pdat_5_5 = melt(pdat[,grep("_5_5|true", colnames(pdat))], id.vars = "true_gen_lia", variable.name = "method")
p = ggplot(pdat_5_5, aes(x = true_gen_lia, y = value, color = method)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(. ~ method) +
  theme_minimal()
plot(p)

betas = true$lia_betas[1:10000] + true$male_betas[1:10000]


files = str_subset(all.files, "v1_")
files = str_subset(files, "_5_5")
files = str_subset(files, "qassoc|logistic")
sig.thres = 5e-8

beta.names = sapply(str_split(files, pattern = "\\_10male_|\\."), FUN = function(x) x[2])
beta.est = data.frame("true.betas" = betas)
p.est = data.frame("true.betas" = betas)
for (i in 1:length(files)) {
  betas.ph = fread(files[i])
  if ("BETA" %in% colnames(betas.ph)) {
    beta.est[[beta.names[i]]] = betas.ph$BETA 
    #  beta.est[[paste(beta.names[i],"_col", sep = "")]] = c("black","red")[(betas.ph$P <= sig.thres) + 1]
  } else {
    beta.est[[beta.names[i]]] = log(betas.ph$OR)
    # beta.est[[paste(beta.names[i],"_col", sep = "")]] = c("black","red")[(betas.ph$P <= sig.thres) + 1]
  }
}
beta.est = melt(beta.est, id.vars = "true.betas", measure.vars = grep(pattern = "_5_5",colnames(beta.est)),
                variable.name = "method")
beta.est[["sig"]] = (betas.ph$P <= sig.thres) + 1


p = ggplot(beta.est, aes(x = true.betas, y = value,color = sig)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(. ~ method) +
  theme_minimal()
plot(p)





#














# Beware - retired code ahead ! -------------------------------------------



library(viridis)
library(hrbrthemes)
p2 = ggplot(bs.res, aes(x = rownames(bs.res), y = mean_C_snps, fill = offset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_C_snps - mean_C_snps_sd, ymax = mean_C_snps + mean_C_snps_sd)) +
  facet_wrap(~type) +
  xlab("") +
  
  theme_minimal()
plot(p2)



imp.vals = lapply(val.list, FUN = function(x) x[,1])
imp.vals = matrix(unlist(imp.vals), ncol = 4, byrow = T)
imp.vals = imp.vals[,2] / imp.vals[,1] - 1

bs.imp = replicate(100000, sample(imp.vals, replace = T))
bs.means.imp = colMeans(bs.imp)
imp.vals.null.means = mean(bs.means.imp)
imp.vals.null.sd = sd(bs.means.imp)



c_ltfh = ltfhsib$CHISQ_LINREG[datrue$lia.betas != 0]
bs = replicate(10000, mean(sample(c_ltfh, replace = T) > 15))
mean(bs) + c(-1,0,1)*sd(bs) * 1.96
hist(bs)

#number of sigg snps
sum(as.numeric(da$P_LINREG) <= 5e-08)


#no sigg snps:
sum(ltfhsib$P_LINREG <= 5e-08)




#ltfh sum stats:
newltfh = fread("D:\\data\\workdata\\707293\\FIUN\\SimData\\bedFiles\\test_10k2_C50_v1_new_ltfh_nosib_00os")

#no sigg snps:
sum(newltfh$P_LINREG <= 5e-08)
#prop of true sigg snps:
(sumnewltfh = sum((newltfh$P_LINREG <= 5e-08) & (datrue$lia_betas != 0)))
#mean of true non sigg snps
mean(newltfh$CHISQ_LINREG[datrue$lia_betas == 0])
#mean of true sigg snps
mean(newltfh$CHISQ_LINREG[datrue$lia_betas != 0])

diffs = which((newltfh$P_LINREG >= 5e-08) & (ltfhsib$P_LINREG <= 5e-08))
datrue[diffs,]
ltfhsib[diffs,]
newltfh[diffs,]
diffs = which((newltfh$P_LINREG < 5e-08) & (ltfhsib$P_LINREG >= 5e-08))
datrue[diffs,]
ltfhsib[diffs,]
newltfh[diffs,]

datue = fread(true.files[1])
newltfh = fread(ltfhNew.files[1])
gwas = fread(gwas.files[1])
gwax = fread(gwax.files[1])
ltfh = fread(ltfh.files[1])

c_da = da[datrue$lia_betas != 0,]
c_ltfhsib = ltfhsib[datrue$lia_betas != 0,]
c_newltfh = newltfh[datrue$lia_betas != 0,]

c_pval = data.frame("GWAS" = gwas[datrue$lia_betas != 0,10], "GWAX" = gwax[datrue$lia_betas != 0,10] , "ltfh" = ltfh[datrue$lia_betas != 0, 10], "new_ltfh" =  newltfh[datrue$lia_betas != 0,10])
colnames(c_pval) = c("GWAS" , "GWAX","ltfh","newltfh")

cor(c_da$P_LINREG, c_newltfh$P_LINREG)
cor(c_ltfhsib$P_LINREG, c_da$P_LINREG)
cor(c_ltfhsib$P_LINREG, c_newltfh$P_LINREG)


par(mfrow = c(2,3))
#pvals plotted against each other
for (i in 1:(dim(c_pval)[2] - 1)) {
  for (j in (i + 1):dim(c_pval)[2]) {
    
    plot(-log10(c_pval[,i]), -log10(c_pval[,j]), pch = 19, xlim = c( 0 ,-log10(5e-08) + 10), ylim = c( 0 ,-log10(5e-08) + 10),
         xlab = colnames(c_pval)[i], ylab = colnames(c_pval)[j], main = "-log10(p)")
    points(-log10(c_pval[,i])[(-log10(c_pval[,i]) >= -log10(5e-08) ) & (-log10(c_pval[,j]) < -log10(5e-08)) ],
           -log10(c_pval[,j])[(-log10(c_pval[,i]) >= -log10(5e-08) ) & (-log10(c_pval[,j]) < -log10(5e-08)) ]
           , pch = 19, col = "red")
    points(-log10(c_pval[,i])[(-log10(c_pval[,i]) < -log10(5e-08) ) & (-log10(c_pval[,j]) >= -log10(5e-08)) ],
           -log10(c_pval[,j])[(-log10(c_pval[,i]) < -log10(5e-08) ) & (-log10(c_pval[,j]) >= -log10(5e-08)) ]
           , pch = 19, col = "green")
    abline(b = 1, a = 0, col = "red")
    abline(h = -log10(5e-08), lty = "dashed")
    abline(v = -log10(5e-08), lty = "dashed")
    
  }
}
par(mfrow = c(1,1))

sumltfh/500
sumnewltfh/500
sumcc/500
(sumltfh - sumcc)/sumcc


c_datrue = datrue[datrue$lia.betas != 0, 1]


#find a way to plot the estimated effect sizes aginst the true effect sizes.
which(c_ltfhsib$P_LINREG <= 5e-8)
cols = rep("black", 500)
cols[which(c_ltfhsib$P_LINREG <= 5e-8)] = "red"
plot(abs(c_datrue$lia.betas), abs(qnorm(c_ltfhsib$P_LINREG) / sqrt(100000)), pch = 19, col = cols,
     xlab = "true betas", ylab = "simulated betas", main = "ltfh_nosib")
abline(h = abs(qnorm(5e-8)/ sqrt(100000)))
mean(c_datrue$lia.betas)
var(c_datrue$lia.betas)
sd(c_datrue$lia.betas)


which(c_da$P_LINREG <= 5e-8)
cols = rep("black", 500)
cols[which(c_da$P_LINREG <= 5e-8)] = "red"
plot(abs(c_datrue$lia.betas), abs(qnorm(c_da$P_LINREG) / sqrt(100000)), pch = 19, col = cols,
     xlab = "true betas", ylab = "simulated betas", main = "GWAS")
abline(h = abs(qnorm(5e-8)/ sqrt(100000)))

par(mfrow = c(2,2))
for (name in colnames(c_pval)) {
  cols = rep("black", 500)
  cols[which(c_pval[[name]] <= 5e-8)] = "red"
  plot(abs(c_datrue$lia.betas), abs(qnorm(c_pval[[name]]) / sqrt(100000)), pch = 19, col = cols,
       xlab = "true betas", ylab = "simulated betas", main = name)
  abline(h = abs(qnorm(5e-8) / sqrt(100000)))
}
par(mfrow = c(1,1))

#####




test = replicate(10000, rchisq(500, df = 33))

test2 = apply(test, MARGIN = 1, FUN = function(x) -log10(pchisq(x, df = 1, lower.tail = F)))
test[,1]
test2[,1]
mean(test[,1])
mean(test2[,1])

test3 = apply(test2, MARGIN = 1, FUN = function(x) x >= 5e-6)
mean(colSums(test3))


hist(ltfhsib$CHISQ_LINREG[datrue$lia.betas != 0], breaks = 100, freq = F)
curve(dchisq(x, df = 33), from = 0, to = 200, n = 1000, add = T, col = "red")
abline(v = cutoff, col = "blue", lty = "dashed")
sum(ltfhsib$CHISQ_LINREG[datrue$lia.betas != 0] > cutoff)
mean(ltfhsib$CHISQ_LINREG[datrue$lia.betas != 0])


datrue[which.max(ltfhsib$CHISQ_LINREG),]

datrue = fread("D:\\data\\workdata\\707293\\FIUN\\SimData\\pedFiles\\noEnvCorrOffPs10kx10k_v1.true")
#ltfh sum stats:
pltfh = fread("D:\\data\\workdata\\707293\\FIUN\\SimData\\bedFiles\\noEnvCorrOffPs10kx10k_v1_ltfh_nosib")
#no sigg snps:
sum(pltfh$P <= 5e-08)
#prop of true sigg snps:
(sumpltfh = sum((pltfh$P <= 5e-08) & (datrue$lia.betas != 0)))
#mean of true non sigg snps
mean(pltfh$CHISQ_LINREG[datrue$lia.betas == 0])
#mean of true sigg snps
mean(pltfh$CHISQ_LINREG[datrue$lia.betas != 0])


#ltfh sum stats:
pda = fread("D:\\data\\workdata\\707293\\FIUN\\SimData\\bedFiles\\noEnvCorrOffPs10kx10k_v1_CHILD_STATUS")
#no sigg snps:
sum(pda$P <= 5e-08)
#prop of true sigg snps:
(sumpda = sum((pda$P <= 5e-08) & (datrue$lia.betas != 0)))
#mean of true non sigg snps
mean(pda$CHISQ[datrue$lia.betas == 0])
#mean of true sigg snps
mean(pda$CHISQ[datrue$lia.betas != 0])







#names(datrue)
#for (i in 1:length(true.files)) {
#  datrue = fread(true.files[i])
##  names(datrue) = gsub("\\_\\_", "\\_", gsub("\\.", "\\_", names(datrue)))
##  names(datrue) = gsub("_cc", "_grp", names(datrue))
#  datrue = cbind("FID" = 1:10000, "IID" = 1:10000, datrue)
#  fwrite(datrue, file = true.files[i], quote = F, sep = " ")
#}
#da1 = fread("D:\\data\\workdata\\707293\\FIUN\\SimData\\pedFiles\\LTFHtest_10k2_C50_v1.phen")
#da2 = fread("D:\\data\\workdata\\707293\\FIUN\\SimData\\pedFiles\\LTFH00ostest_10k2_C50_v1.phen")
#table(da1$ltfh, da2$ltfh)
